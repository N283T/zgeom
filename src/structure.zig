const std = @import("std");
const geometry = @import("geometry.zig");

pub const Atom = struct {
    model: u32,
    chain: []const u8,
    seq: []const u8,
    ins: []const u8,
    residue: []const u8,
    name: []const u8,
    alt: []const u8,
    occupancy: f64,
    pos: geometry.Vec3,
    polymer_hint: bool = false,
    segment: u32 = 0,
};

pub const Structure = struct {
    allocator: std.mem.Allocator,
    source: []u8,
    atoms: []Atom,

    pub fn deinit(self: *Structure) void {
        self.allocator.free(self.atoms);
        self.allocator.free(self.source);
        self.* = undefined;
    }

    pub fn firstModel(self: Structure) ?u32 {
        if (self.atoms.len == 0) return null;
        var result = self.atoms[0].model;
        for (self.atoms[1..]) |atom| result = @min(result, atom.model);
        return result;
    }
};

pub const AtomSpec = struct {
    chain: []const u8,
    seq: []const u8,
    ins: []const u8,
    name: []const u8,
    alt: ?[]const u8,

    pub fn parse(text: []const u8) !AtomSpec {
        var parts = std.mem.splitScalar(u8, text, ':');
        const chain = parts.next() orelse return error.InvalidAtomSpec;
        const seq_ins = parts.next() orelse return error.InvalidAtomSpec;
        const atom_alt = parts.next() orelse return error.InvalidAtomSpec;
        if (parts.next() != null or seq_ins.len == 0 or atom_alt.len == 0) return error.InvalidAtomSpec;

        var ins: []const u8 = "";
        var seq = seq_ins;
        if (std.mem.indexOfScalar(u8, seq_ins, '^')) |idx| {
            seq = seq_ins[0..idx];
            ins = seq_ins[idx + 1 ..];
            if (seq.len == 0 or ins.len == 0) return error.InvalidAtomSpec;
        }
        var alt: ?[]const u8 = null;
        var name = atom_alt;
        if (std.mem.indexOfScalar(u8, atom_alt, '@')) |idx| {
            name = atom_alt[0..idx];
            const value = atom_alt[idx + 1 ..];
            if (name.len == 0 or value.len == 0) return error.InvalidAtomSpec;
            alt = value;
        }
        return .{ .chain = chain, .seq = seq, .ins = ins, .name = name, .alt = alt };
    }
};

fn slice(line: []const u8, start: usize, end: usize) []const u8 {
    if (start >= line.len) return "";
    return line[start..@min(end, line.len)];
}

fn cleanCif(value: []const u8) []const u8 {
    if (std.mem.eql(u8, value, ".") or std.mem.eql(u8, value, "?")) return "";
    return value;
}

fn cifEql(a: []const u8, b: []const u8) bool {
    return std.ascii.eqlIgnoreCase(a, b);
}

fn cifStartsWith(value: []const u8, prefix: []const u8) bool {
    return value.len >= prefix.len and std.ascii.eqlIgnoreCase(value[0..prefix.len], prefix);
}

fn isRecognizedProteinResidue(name: []const u8) bool {
    const names = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "UNK", "MSE", "SEC", "PYL",
    };
    for (names) |candidate| if (std.mem.eql(u8, name, candidate)) return true;
    return false;
}

fn parseFiniteFloat(value: []const u8) !f64 {
    const result = std.fmt.parseFloat(f64, value) catch return error.InvalidStructureNumber;
    if (!std.math.isFinite(result)) return error.InvalidStructureNumber;
    return result;
}

fn stableId(value: []const u8) u32 {
    var result: u32 = 2166136261;
    for (value) |byte| result = (result ^ byte) *% 16777619;
    return result;
}

pub fn load(io: std.Io, allocator: std.mem.Allocator, path: []const u8) !Structure {
    const raw = try std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(1024 * 1024 * 1024));
    const source = if (raw.len >= 2 and raw[0] == 0x1f and raw[1] == 0x8b) decompress: {
        errdefer allocator.free(raw);
        var source_reader: std.Io.Reader = .fixed(raw);
        var window: [std.compress.flate.max_window_len]u8 = undefined;
        var decompressor: std.compress.flate.Decompress = .init(&source_reader, .gzip, &window);
        var output: std.Io.Writer.Allocating = .init(allocator);
        errdefer output.deinit();
        _ = decompressor.reader.streamRemaining(&output.writer) catch return error.InvalidCompressedInput;
        if (output.writer.buffered().len > 1024 * 1024 * 1024) return error.StreamTooLong;
        const decompressed = try output.toOwnedSlice();
        allocator.free(raw);
        break :decompress decompressed;
    } else raw;
    errdefer allocator.free(source);
    var atoms = std.ArrayListUnmanaged(Atom).empty;
    errdefer atoms.deinit(allocator);

    const trimmed = std.mem.trimStart(u8, source, " \t\r\n");
    if (std.mem.endsWith(u8, path, ".pdb") or std.mem.endsWith(u8, path, ".pdb.gz") or std.mem.startsWith(u8, trimmed, "ATOM  ") or std.mem.startsWith(u8, trimmed, "HEADER")) {
        try parsePdb(allocator, source, &atoms);
    } else if (std.mem.endsWith(u8, path, ".cif") or std.mem.endsWith(u8, path, ".cif.gz") or std.mem.startsWith(u8, trimmed, "data_")) {
        try parseMmcif(allocator, source, &atoms);
    } else return error.UnsupportedFormat;
    if (atoms.items.len == 0) return error.NoAtoms;
    return .{ .allocator = allocator, .source = source, .atoms = try atoms.toOwnedSlice(allocator) };
}

fn parsePdb(allocator: std.mem.Allocator, source: []const u8, atoms: *std.ArrayListUnmanaged(Atom)) !void {
    var current_model: u32 = 1;
    var segment: u32 = 0;
    var lines = std.mem.splitScalar(u8, source, '\n');
    while (lines.next()) |raw_line| {
        const line = std.mem.trimEnd(u8, raw_line, "\r");
        if (std.mem.startsWith(u8, line, "MODEL ")) {
            current_model = std.fmt.parseInt(u32, std.mem.trim(u8, slice(line, 10, 14), " "), 10) catch return error.InvalidModel;
            continue;
        }
        if (std.mem.startsWith(u8, line, "TER")) {
            segment +%= 1;
            continue;
        }
        if (!std.mem.startsWith(u8, line, "ATOM  ") and !std.mem.startsWith(u8, line, "HETATM")) continue;
        if (line.len < 54) return error.TruncatedPdbAtom;
        const x = try parseFiniteFloat(std.mem.trim(u8, slice(line, 30, 38), " "));
        const y = try parseFiniteFloat(std.mem.trim(u8, slice(line, 38, 46), " "));
        const z = try parseFiniteFloat(std.mem.trim(u8, slice(line, 46, 54), " "));
        const occupancy_text = std.mem.trim(u8, slice(line, 54, 60), " ");
        const occupancy = if (occupancy_text.len == 0) 1.0 else try parseFiniteFloat(occupancy_text);
        const residue_name = std.mem.trim(u8, slice(line, 17, 20), " ");
        try atoms.append(allocator, .{
            .model = current_model,
            .chain = std.mem.trim(u8, slice(line, 21, 22), " "),
            .seq = std.mem.trim(u8, slice(line, 22, 26), " "),
            .ins = std.mem.trim(u8, slice(line, 26, 27), " "),
            .residue = residue_name,
            .name = std.mem.trim(u8, slice(line, 12, 16), " "),
            .alt = std.mem.trim(u8, slice(line, 16, 17), " "),
            .occupancy = occupancy,
            .pos = .{ .x = x, .y = y, .z = z },
            .polymer_hint = isRecognizedProteinResidue(residue_name),
            .segment = segment,
        });
    }
}

fn tokenizeCif(allocator: std.mem.Allocator, source: []const u8) ![][]const u8 {
    var tokens = std.ArrayListUnmanaged([]const u8).empty;
    errdefer tokens.deinit(allocator);
    var i: usize = 0;
    var line_start = true;
    while (i < source.len) {
        if (source[i] == '#') {
            while (i < source.len and source[i] != '\n') i += 1;
            line_start = true;
            continue;
        }
        if (std.ascii.isWhitespace(source[i])) {
            line_start = source[i] == '\n';
            i += 1;
            continue;
        }
        if (source[i] == ';' and line_start) {
            const start = i + 1;
            i = start;
            while (i + 1 < source.len and !(source[i] == '\n' and source[i + 1] == ';')) i += 1;
            try tokens.append(allocator, source[start..i]);
            if (i < source.len) {
                i += 2;
                while (i < source.len and source[i] != '\n') i += 1;
            }
            line_start = true;
            continue;
        }
        line_start = false;
        if (source[i] == '\'' or source[i] == '"') {
            const quote = source[i];
            i += 1;
            const start = i;
            while (i < source.len) : (i += 1) {
                if (source[i] == quote and (i + 1 == source.len or std.ascii.isWhitespace(source[i + 1]))) break;
            }
            if (i >= source.len) return error.UnterminatedCifQuote;
            try tokens.append(allocator, source[start..i]);
            i += 1;
        } else {
            const start = i;
            while (i < source.len and !std.ascii.isWhitespace(source[i])) : (i += 1) {}
            try tokens.append(allocator, source[start..i]);
        }
    }
    return tokens.toOwnedSlice(allocator);
}

fn findTag(tags: []const []const u8, suffixes: []const []const u8) ?usize {
    for (suffixes) |suffix| for (tags, 0..) |tag, idx| if (cifEql(tag, suffix)) return idx;
    return null;
}

fn parseMmcif(allocator: std.mem.Allocator, source: []const u8, atoms: *std.ArrayListUnmanaged(Atom)) !void {
    const tokens = try tokenizeCif(allocator, source);
    defer allocator.free(tokens);
    var i: usize = 0;
    while (i < tokens.len) : (i += 1) {
        if (!cifEql(tokens[i], "loop_")) continue;
        i += 1;
        const tag_start = i;
        while (i < tokens.len and std.mem.startsWith(u8, tokens[i], "_")) i += 1;
        const tags = tokens[tag_start..i];
        if (tags.len == 0 or !cifStartsWith(tags[0], "_atom_site.")) {
            i -= 1;
            continue;
        }
        const group_i = findTag(tags, &.{"_atom_site.group_PDB"});
        const atom_i = findTag(tags, &.{ "_atom_site.auth_atom_id", "_atom_site.label_atom_id" }) orelse return error.MissingCifColumn;
        const alt_i = findTag(tags, &.{"_atom_site.label_alt_id"});
        const comp_i = findTag(tags, &.{ "_atom_site.auth_comp_id", "_atom_site.label_comp_id" }) orelse return error.MissingCifColumn;
        const chain_i = findTag(tags, &.{ "_atom_site.auth_asym_id", "_atom_site.label_asym_id" }) orelse return error.MissingCifColumn;
        const label_chain_i = findTag(tags, &.{"_atom_site.label_asym_id"});
        const seq_i = findTag(tags, &.{ "_atom_site.auth_seq_id", "_atom_site.label_seq_id" }) orelse return error.MissingCifColumn;
        const label_seq_i = findTag(tags, &.{"_atom_site.label_seq_id"});
        const ins_i = findTag(tags, &.{"_atom_site.pdbx_PDB_ins_code"});
        const x_i = findTag(tags, &.{"_atom_site.Cartn_x"}) orelse return error.MissingCifColumn;
        const y_i = findTag(tags, &.{"_atom_site.Cartn_y"}) orelse return error.MissingCifColumn;
        const z_i = findTag(tags, &.{"_atom_site.Cartn_z"}) orelse return error.MissingCifColumn;
        const occ_i = findTag(tags, &.{"_atom_site.occupancy"});
        const model_i = findTag(tags, &.{"_atom_site.pdbx_PDB_model_num"});
        while (i + tags.len <= tokens.len) : (i += tags.len) {
            if (cifEql(tokens[i], "loop_") or std.mem.startsWith(u8, tokens[i], "_") or cifStartsWith(tokens[i], "data_") or cifEql(tokens[i], "stop_")) break;
            const row = tokens[i .. i + tags.len];
            if (group_i) |idx| if (!cifEql(row[idx], "ATOM") and !cifEql(row[idx], "HETATM")) continue;
            try atoms.append(allocator, .{
                .model = if (model_i) |idx| try std.fmt.parseInt(u32, row[idx], 10) else 1,
                .chain = cleanCif(row[chain_i]),
                .seq = cleanCif(row[seq_i]),
                .ins = if (ins_i) |idx| cleanCif(row[idx]) else "",
                .residue = cleanCif(row[comp_i]),
                .name = cleanCif(row[atom_i]),
                .alt = if (alt_i) |idx| cleanCif(row[idx]) else "",
                .occupancy = if (occ_i) |idx| if (cleanCif(row[idx]).len == 0) 1.0 else try parseFiniteFloat(row[idx]) else 1.0,
                .pos = .{
                    .x = try parseFiniteFloat(row[x_i]),
                    .y = try parseFiniteFloat(row[y_i]),
                    .z = try parseFiniteFloat(row[z_i]),
                },
                .polymer_hint = (if (label_seq_i) |idx| cleanCif(row[idx]).len != 0 else true) and isRecognizedProteinResidue(cleanCif(row[comp_i])),
                .segment = if (label_chain_i) |idx| stableId(cleanCif(row[idx])) else stableId(cleanCif(row[chain_i])),
            });
        }
        if (i < tokens.len and
            !cifEql(tokens[i], "loop_") and
            !std.mem.startsWith(u8, tokens[i], "_") and
            !cifStartsWith(tokens[i], "data_") and
            !cifEql(tokens[i], "stop_")) return error.IncompleteCifRow;
        if (i > 0) i -= 1;
    }
}

fn sameIdentity(atom: Atom, model: u32, chain: []const u8, seq: []const u8, ins: []const u8, name: []const u8) bool {
    return atom.model == model and std.mem.eql(u8, atom.chain, chain) and std.mem.eql(u8, atom.seq, seq) and std.mem.eql(u8, atom.ins, ins) and std.mem.eql(u8, atom.name, name);
}

fn betterAlt(candidate: Atom, current: Atom) bool {
    if (candidate.alt.len == 0 and current.alt.len != 0) return true;
    if (candidate.alt.len != 0 and current.alt.len == 0) return false;
    if (candidate.occupancy != current.occupancy) return candidate.occupancy > current.occupancy;
    if (std.mem.eql(u8, candidate.alt, "A") and !std.mem.eql(u8, current.alt, "A")) return true;
    return std.mem.order(u8, candidate.alt, current.alt) == .lt;
}

pub fn findAtom(structure: Structure, model: u32, spec: AtomSpec, global_alt: ?[]const u8) !*const Atom {
    var best: ?*const Atom = null;
    for (structure.atoms) |*atom| {
        if (!sameIdentity(atom.*, model, spec.chain, spec.seq, spec.ins, spec.name)) continue;
        if (spec.alt) |alt| {
            if (!std.mem.eql(u8, atom.alt, alt)) continue;
        } else if (global_alt) |alt| {
            // A requested conformer is combined with atoms that have no
            // alternate identifier, as required to obtain a complete model.
            if (atom.alt.len != 0 and !std.mem.eql(u8, atom.alt, alt)) continue;
        }
        if (best == null or betterAlt(atom.*, best.?.*)) best = atom;
    }
    return best orelse error.AtomNotFound;
}

pub const Residue = struct {
    model: u32,
    chain: []const u8,
    seq: []const u8,
    ins: []const u8,
    name: []const u8,
    atom_start: usize,
    atom_end: usize,
    segment: u32,
};

const ResidueKey = struct {
    segment: u32,
    chain: []const u8,
    seq: []const u8,
    ins: []const u8,
};

const ResidueKeyContext = struct {
    pub fn hash(_: @This(), key: ResidueKey) u64 {
        var hasher = std.hash.Wyhash.init(key.segment);
        hasher.update(key.chain);
        hasher.update(&.{0});
        hasher.update(key.seq);
        hasher.update(&.{0});
        hasher.update(key.ins);
        return hasher.final();
    }

    pub fn eql(_: @This(), a: ResidueKey, b: ResidueKey) bool {
        return a.segment == b.segment and
            std.mem.eql(u8, a.chain, b.chain) and
            std.mem.eql(u8, a.seq, b.seq) and
            std.mem.eql(u8, a.ins, b.ins);
    }
};

fn residueKey(atom: Atom) ResidueKey {
    return .{ .segment = atom.segment, .chain = atom.chain, .seq = atom.seq, .ins = atom.ins };
}

fn validateContiguousResidues(structure: Structure, allocator: std.mem.Allocator, model: u32) !void {
    var seen: std.HashMapUnmanaged(ResidueKey, void, ResidueKeyContext, 80) = .empty;
    defer seen.deinit(allocator);
    var current: ?ResidueKey = null;
    for (structure.atoms) |atom| {
        if (atom.model != model or !atom.polymer_hint) continue;
        const key = residueKey(atom);
        if (current) |previous| {
            if (ResidueKeyContext.eql(.{}, previous, key)) continue;
        }
        const entry = try seen.getOrPut(allocator, key);
        if (entry.found_existing) return error.NonContiguousResidue;
        current = key;
    }
}

fn preferredAlt(candidate: []const u8, current: []const u8) bool {
    if (std.mem.eql(u8, candidate, "A") and !std.mem.eql(u8, current, "A")) return true;
    if (!std.mem.eql(u8, candidate, "A") and std.mem.eql(u8, current, "A")) return false;
    return std.mem.order(u8, candidate, current) == .lt;
}

/// Select one named conformer for an entire residue. Blank-altloc atoms are
/// shared by every conformer and therefore do not participate in the score.
/// The default uses the greatest mean occupancy across rows carrying each
/// named altloc; exact ties prefer A, then lexical order. An explicit request
/// is used wherever the residue has named alternates; blank-only residues
/// report a blank label because all of their atoms are shared.
pub fn resolveResidueAlt(structure: Structure, residue: Residue, requested: ?[]const u8) []const u8 {
    var best_alt: []const u8 = "";
    var best_total: f64 = 0.0;
    var best_count: usize = 0;
    const atoms = structure.atoms[residue.atom_start..residue.atom_end];
    if (requested) |alt| {
        // A residue containing only blank/shared atoms has no alternate
        // conformer to report. If named alternatives do exist, retain the
        // request even when that label is absent: lookups then expose missing
        // coordinates rather than silently falling back to another label.
        for (atoms) |atom| if (atom.alt.len != 0) return alt;
        return "";
    }
    for (atoms, 0..) |candidate, candidate_idx| {
        if (candidate.alt.len == 0) continue;
        var already_scored = false;
        for (atoms[0..candidate_idx]) |earlier| {
            if (std.mem.eql(u8, earlier.alt, candidate.alt)) {
                already_scored = true;
                break;
            }
        }
        if (already_scored) continue;

        var total: f64 = 0.0;
        var count: usize = 0;
        for (atoms) |atom| {
            if (!std.mem.eql(u8, atom.alt, candidate.alt)) continue;
            total += atom.occupancy;
            count += 1;
        }
        if (best_count == 0) {
            best_alt = candidate.alt;
            best_total = total;
            best_count = count;
            continue;
        }
        const candidate_weighted = total * @as(f64, @floatFromInt(best_count));
        const best_weighted = best_total * @as(f64, @floatFromInt(count));
        if (candidate_weighted > best_weighted or
            (candidate_weighted == best_weighted and preferredAlt(candidate.alt, best_alt)))
        {
            best_alt = candidate.alt;
            best_total = total;
            best_count = count;
        }
    }
    return best_alt;
}

/// Resolve one atom inside a residue discovered by `residues`. Restricting the
/// scan to the residue's contiguous atom range keeps backbone analysis linear
/// in structure size rather than rescanning all atoms for every torsion.
pub fn findResidueAtom(structure: Structure, residue: Residue, name: []const u8, global_alt: ?[]const u8) !*const Atom {
    var best: ?*const Atom = null;
    for (structure.atoms[residue.atom_start..residue.atom_end]) |*atom| {
        if (atom.model != residue.model or
            !std.mem.eql(u8, atom.chain, residue.chain) or
            !std.mem.eql(u8, atom.seq, residue.seq) or
            !std.mem.eql(u8, atom.ins, residue.ins) or
            atom.segment != residue.segment or
            !std.mem.eql(u8, atom.name, name)) continue;
        if (global_alt) |alt| {
            if (atom.alt.len != 0 and !std.mem.eql(u8, atom.alt, alt)) continue;
        }
        if (best == null or betterAlt(atom.*, best.?.*)) best = atom;
    }
    return best orelse error.AtomNotFound;
}

pub fn residues(structure: Structure, allocator: std.mem.Allocator, model: u32, chain_filter: ?[]const u8) ![]Residue {
    try validateContiguousResidues(structure, allocator, model);
    var result = std.ArrayListUnmanaged(Residue).empty;
    errdefer result.deinit(allocator);
    for (structure.atoms, 0..) |atom, atom_idx| {
        if (atom.model != model) continue;
        if (!atom.polymer_hint) continue;
        if (chain_filter) |chain| if (!std.mem.eql(u8, atom.chain, chain)) continue;
        if (result.items.len > 0) {
            const previous = &result.items[result.items.len - 1];
            if (previous.segment == atom.segment and std.mem.eql(u8, previous.chain, atom.chain) and std.mem.eql(u8, previous.seq, atom.seq) and std.mem.eql(u8, previous.ins, atom.ins)) {
                previous.atom_end = atom_idx + 1;
                continue;
            }
        }
        try result.append(allocator, .{
            .model = model,
            .chain = atom.chain,
            .seq = atom.seq,
            .ins = atom.ins,
            .name = atom.residue,
            .atom_start = atom_idx,
            .atom_end = atom_idx + 1,
            .segment = atom.segment,
        });
    }
    return result.toOwnedSlice(allocator);
}

test "PDB parsing, models, insertion code and altloc policy" {
    const allocator = std.testing.allocator;
    const text =
        "MODEL        2\n" ++
        "ATOM      1  CA AALA A  10A      1.000   2.000   3.000  0.40 20.00           C\n" ++
        "ATOM      2  CA BALA A  10A      4.000   5.000   6.000  0.60 20.00           C\n" ++
        "ENDMDL\n";
    const source = try allocator.dupe(u8, text);
    var list = std.ArrayListUnmanaged(Atom).empty;
    try parsePdb(allocator, source, &list);
    var structure = Structure{ .allocator = allocator, .source = source, .atoms = try list.toOwnedSlice(allocator) };
    defer structure.deinit();
    const spec = try AtomSpec.parse("A:10^A:CA");
    const atom = try findAtom(structure, 2, spec, null);
    try std.testing.expectEqualStrings("B", atom.alt);
    const atom_a = try findAtom(structure, 2, spec, "A");
    try std.testing.expectEqualStrings("A", atom_a.alt);
}

test "mmCIF atom_site parsing" {
    const allocator = std.testing.allocator;
    const text = "data_x\nloop_\n_atom_site.group_PDB\n_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.pdbx_PDB_model_num\nATOM CA . GLY A 1 ? 1.0 2.0 3.0 1.0 1\n#\n";
    var list = std.ArrayListUnmanaged(Atom).empty;
    defer list.deinit(allocator);
    try parseMmcif(allocator, text, &list);
    try std.testing.expectEqual(@as(usize, 1), list.items.len);
    try std.testing.expectEqualStrings("CA", list.items[0].name);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), list.items[0].pos.y, 1e-12);
}

test "load stable PDB and mmCIF fixtures" {
    const allocator = std.testing.allocator;
    var pdb = try load(std.testing.io, allocator, "tests/fixtures/edge_cases.pdb");
    defer pdb.deinit();
    try std.testing.expectEqual(@as(usize, 10), pdb.atoms.len);
    try std.testing.expectEqual(@as(?u32, 1), pdb.firstModel());

    var cif = try load(std.testing.io, allocator, "tests/fixtures/simple.cif");
    defer cif.deinit();
    try std.testing.expectEqual(@as(usize, 7), cif.atoms.len);
    const atom = try findAtom(cif, 1, try AtomSpec.parse("A:2:CA"), null);
    try std.testing.expectApproxEqAbs(@as(f64, 2.3), atom.pos.x, 1e-12);
    const residue_list = try residues(cif, allocator, 1, null);
    defer allocator.free(residue_list);
    try std.testing.expectEqual(@as(usize, 2), residue_list.len);
}

test "residue-level altloc selection avoids mixed backbone conformers" {
    const allocator = std.testing.allocator;
    var pdb = try load(std.testing.io, allocator, "tests/fixtures/coherent_altloc.pdb");
    defer pdb.deinit();
    const residue_list = try residues(pdb, allocator, 1, null);
    defer allocator.free(residue_list);

    try std.testing.expectEqual(@as(usize, 2), residue_list.len);
    const first_alt = resolveResidueAlt(pdb, residue_list[0], null);
    try std.testing.expectEqualStrings("A", first_alt);
    try std.testing.expectEqualStrings("A", (try findResidueAtom(pdb, residue_list[0], "N", first_alt)).alt);
    try std.testing.expectEqualStrings("A", (try findResidueAtom(pdb, residue_list[0], "CA", first_alt)).alt);
    try std.testing.expectEqualStrings("A", (try findResidueAtom(pdb, residue_list[0], "C", first_alt)).alt);

    try std.testing.expectEqualStrings("B", resolveResidueAlt(pdb, residue_list[0], "B"));
    try std.testing.expectEqualStrings("", resolveResidueAlt(pdb, residue_list[1], "B"));
}

test "non-contiguous residue atom rows are rejected" {
    const allocator = std.testing.allocator;
    var cif = try load(std.testing.io, allocator, "tests/fixtures/noncontiguous_residue.cif");
    defer cif.deinit();
    try std.testing.expectError(error.NonContiguousResidue, residues(cif, allocator, 1, null));
}

test "recognized modified amino acid HETATM rows remain protein backbone" {
    const allocator = std.testing.allocator;
    var pdb = try load(std.testing.io, allocator, "tests/fixtures/modified_protein.pdb");
    defer pdb.deinit();
    const residue_list = try residues(pdb, allocator, 1, null);
    defer allocator.free(residue_list);
    try std.testing.expectEqual(@as(usize, 2), residue_list.len);
    try std.testing.expectEqualStrings("MSE", residue_list[0].name);
    try std.testing.expectEqualStrings("GLY", residue_list[1].name);
}

test "corrupt gzip input is a malformed structure error" {
    try std.testing.expectError(
        error.InvalidCompressedInput,
        load(std.testing.io, std.testing.allocator, "tests/fixtures/corrupt.cif.gz"),
    );
}

test "blank chain atom spec and explicit PDB TER segments" {
    const blank = try AtomSpec.parse(":1:CA");
    try std.testing.expectEqualStrings("", blank.chain);

    const allocator = std.testing.allocator;
    var pdb = try load(std.testing.io, allocator, "tests/fixtures/ter_break.pdb");
    defer pdb.deinit();
    const residue_list = try residues(pdb, allocator, 1, null);
    defer allocator.free(residue_list);
    try std.testing.expectEqual(@as(usize, 2), residue_list.len);
    try std.testing.expect(residue_list[0].segment != residue_list[1].segment);
}
