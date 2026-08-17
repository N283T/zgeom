const std = @import("std");
const zgeom = @import("zgeom");
const build_options = @import("build_options");

const Format = enum { tsv, csv, json };
const Unit = enum { angstrom, degree, radian };

const Options = struct {
    model: ?u32 = null,
    format: Format = .tsv,
    unit: ?Unit = null,
    chain: ?[]const u8 = null,
    altloc: ?[]const u8 = null,
};

fn usage() void {
    std.debug.print(
        \\zgeom {s} — standalone molecular geometry CLI
        \\
        \\Usage:
        \\  zgeom distance  FILE ATOM1 ATOM2 [OPTIONS]
        \\  zgeom angle     FILE ATOM1 ATOM2 ATOM3 [OPTIONS]
        \\  zgeom dihedral  FILE ATOM1 ATOM2 ATOM3 ATOM4 [OPTIONS]
        \\  zgeom backbone  FILE [OPTIONS]
        \\
        \\Atom syntax: CHAIN:AUTH_SEQ[^INS_CODE]:ATOM_NAME[@ALTLOC]
        \\Examples: A:42:CA  A:42^B:N  A:42:CB@B
        \\
        \\Options:
        \\  --model N              model number (default: lowest model)
        \\  --altloc ID            preferred altloc; blank atoms remain eligible
        \\  --chain ID             restrict backbone output to one chain
        \\  --unit angstrom|deg|rad distance is always angstrom; angles default deg
        \\  --format tsv|csv|json  output format (default: tsv)
        \\  -h, --help             show this help
        \\  -V, --version          show version
        \\
    , .{build_options.version});
}

fn parseOptions(args: []const []const u8, start: usize) !Options {
    var options = Options{};
    var i = start;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "--model")) {
            i += 1;
            if (i >= args.len) return error.MissingOptionValue;
            options.model = std.fmt.parseInt(u32, args[i], 10) catch return error.InvalidOptionValue;
        } else if (std.mem.eql(u8, arg, "--format")) {
            i += 1;
            if (i >= args.len) return error.MissingOptionValue;
            options.format = if (std.mem.eql(u8, args[i], "tsv")) .tsv else if (std.mem.eql(u8, args[i], "csv")) .csv else if (std.mem.eql(u8, args[i], "json")) .json else return error.InvalidOptionValue;
        } else if (std.mem.eql(u8, arg, "--unit")) {
            i += 1;
            if (i >= args.len) return error.MissingOptionValue;
            options.unit = if (std.mem.eql(u8, args[i], "angstrom")) .angstrom else if (std.mem.eql(u8, args[i], "deg")) .degree else if (std.mem.eql(u8, args[i], "rad")) .radian else return error.InvalidOptionValue;
        } else if (std.mem.eql(u8, arg, "--chain")) {
            i += 1;
            if (i >= args.len or args[i].len == 0) return error.MissingOptionValue;
            options.chain = args[i];
        } else if (std.mem.eql(u8, arg, "--altloc")) {
            i += 1;
            if (i >= args.len or args[i].len == 0) return error.MissingOptionValue;
            options.altloc = args[i];
        } else return error.UnknownOption;
    }
    return options;
}

fn resolvedModel(structure: zgeom.structure.Structure, requested: ?u32) !u32 {
    const model = requested orelse structure.firstModel() orelse return error.NoAtoms;
    for (structure.atoms) |atom| if (atom.model == model) return model;
    return error.ModelNotFound;
}

fn writeJsonString(writer: *std.Io.Writer, value: []const u8) !void {
    try writer.writeByte('"');
    for (value) |c| switch (c) {
        '"' => try writer.writeAll("\\\""),
        '\\' => try writer.writeAll("\\\\"),
        '\n' => try writer.writeAll("\\n"),
        '\r' => try writer.writeAll("\\r"),
        '\t' => try writer.writeAll("\\t"),
        else => if (c < 0x20) try writer.print("\\u{x:0>4}", .{c}) else try writer.writeByte(c),
    };
    try writer.writeByte('"');
}

fn unitName(unit: Unit) []const u8 {
    return switch (unit) {
        .angstrom => "angstrom",
        .degree => "degree",
        .radian => "radian",
    };
}

fn outputScalar(writer: *std.Io.Writer, format: Format, kind: []const u8, model: u32, specs: []const []const u8, value: f64, unit: Unit) !void {
    switch (format) {
        .tsv => {
            try writer.writeAll("kind\tmodel");
            for (specs, 0..) |_, idx| try writer.print("\tatom{d}", .{idx + 1});
            try writer.writeAll("\tvalue\tunit\n");
            try writer.print("{s}\t{d}", .{ kind, model });
            for (specs) |spec| try writer.print("\t{s}", .{spec});
            try writer.print("\t{d:.6}\t{s}\n", .{ value, unitName(unit) });
        },
        .csv => {
            try writer.writeAll("kind,model");
            for (specs, 0..) |_, idx| try writer.print(",atom{d}", .{idx + 1});
            try writer.writeAll(",value,unit\n");
            try writer.print("{s},{d}", .{ kind, model });
            for (specs) |spec| try writer.print(",{s}", .{spec});
            try writer.print(",{d:.6},{s}\n", .{ value, unitName(unit) });
        },
        .json => {
            try writer.writeAll("{\"kind\":");
            try writeJsonString(writer, kind);
            try writer.print(",\"model\":{d},\"atoms\":[", .{model});
            for (specs, 0..) |spec, idx| {
                if (idx != 0) try writer.writeByte(',');
                try writeJsonString(writer, spec);
            }
            try writer.print("],\"value\":{d:.6},\"unit\":", .{value});
            try writeJsonString(writer, unitName(unit));
            try writer.writeAll("}\n");
        },
    }
}

fn runScalar(io: std.Io, allocator: std.mem.Allocator, command: []const u8, args: []const []const u8) !void {
    const atom_count: usize = if (std.mem.eql(u8, command, "distance")) 2 else if (std.mem.eql(u8, command, "angle")) 3 else 4;
    if (args.len < 1 + atom_count) return error.MissingArgument;
    const options = try parseOptions(args, 1 + atom_count);
    var structure = try zgeom.structure.load(io, allocator, args[0]);
    defer structure.deinit();
    const model = try resolvedModel(structure, options.model);

    var positions: [4]zgeom.geometry.Vec3 = undefined;
    var labels: [4][]u8 = undefined;
    var labels_initialized: usize = 0;
    defer for (labels[0..labels_initialized]) |label| allocator.free(label);
    for (0..atom_count) |idx| {
        const spec = try zgeom.structure.AtomSpec.parse(args[idx + 1]);
        const atom = try zgeom.structure.findAtom(structure, model, spec, options.altloc);
        positions[idx] = atom.pos;
        labels[idx] = if (atom.ins.len == 0 and atom.alt.len == 0)
            try std.fmt.allocPrint(allocator, "{s}:{s}:{s}", .{ atom.chain, atom.seq, atom.name })
        else if (atom.ins.len == 0)
            try std.fmt.allocPrint(allocator, "{s}:{s}:{s}@{s}", .{ atom.chain, atom.seq, atom.name, atom.alt })
        else if (atom.alt.len == 0)
            try std.fmt.allocPrint(allocator, "{s}:{s}^{s}:{s}", .{ atom.chain, atom.seq, atom.ins, atom.name })
        else
            try std.fmt.allocPrint(allocator, "{s}:{s}^{s}:{s}@{s}", .{ atom.chain, atom.seq, atom.ins, atom.name, atom.alt });
        labels_initialized += 1;
    }
    var unit: Unit = if (atom_count == 2) .angstrom else .degree;
    if (options.unit) |requested| unit = requested;
    if (atom_count == 2 and unit != .angstrom) return error.InvalidUnit;
    if (atom_count > 2 and unit == .angstrom) return error.InvalidUnit;

    const value = if (atom_count == 2)
        zgeom.geometry.distance(positions[0], positions[1])
    else if (atom_count == 3)
        zgeom.geometry.angle(positions[0], positions[1], positions[2]) orelse return error.UndefinedGeometry
    else
        zgeom.geometry.dihedral(positions[0], positions[1], positions[2], positions[3]) orelse return error.UndefinedGeometry;
    const converted = if (unit == .degree) zgeom.geometry.radiansToDegrees(value) else value;

    const stdout = std.Io.File.stdout();
    var buffer: [4096]u8 = undefined;
    var output = stdout.writer(io, &buffer);
    try outputScalar(&output.interface, options.format, command, model, labels[0..atom_count], converted, unit);
    try output.interface.flush();
}

fn residueAtom(structure: zgeom.structure.Structure, res: zgeom.structure.Residue, name: []const u8, altloc: ?[]const u8) ?*const zgeom.structure.Atom {
    return zgeom.structure.findResidueAtom(structure, res, name, altloc) catch null;
}

fn peptideLinked(structure: zgeom.structure.Structure, left: zgeom.structure.Residue, right: zgeom.structure.Residue, left_alt: []const u8, right_alt: []const u8) bool {
    if (left.segment != right.segment or !std.mem.eql(u8, left.chain, right.chain)) return false;
    const c = residueAtom(structure, left, "C", left_alt) orelse return false;
    const n = residueAtom(structure, right, "N", right_alt) orelse return false;
    return zgeom.geometry.distance(c.pos, n.pos) <= 2.0;
}

fn torsion4(a: ?*const zgeom.structure.Atom, b: ?*const zgeom.structure.Atom, c: ?*const zgeom.structure.Atom, d: ?*const zgeom.structure.Atom) ?f64 {
    return zgeom.geometry.dihedral((a orelse return null).pos, (b orelse return null).pos, (c orelse return null).pos, (d orelse return null).pos);
}

fn writeOptional(writer: *std.Io.Writer, value: ?f64, unit: Unit, null_text: []const u8) !void {
    if (value) |v| try writer.print("{d:.6}", .{if (unit == .degree) zgeom.geometry.radiansToDegrees(v) else v}) else try writer.writeAll(null_text);
}

fn runBackbone(io: std.Io, allocator: std.mem.Allocator, args: []const []const u8) !void {
    if (args.len < 1) return error.MissingArgument;
    const options = try parseOptions(args, 1);
    const unit = options.unit orelse .degree;
    if (unit == .angstrom) return error.InvalidUnit;
    var structure = try zgeom.structure.load(io, allocator, args[0]);
    defer structure.deinit();
    const model = try resolvedModel(structure, options.model);
    const residue_list = try zgeom.structure.residues(structure, allocator, model, options.chain);
    defer allocator.free(residue_list);
    if (residue_list.len == 0) return error.NoResidues;
    const selected_alts = try allocator.alloc([]const u8, residue_list.len);
    defer allocator.free(selected_alts);
    for (residue_list, 0..) |res, idx|
        selected_alts[idx] = zgeom.structure.resolveResidueAlt(structure, res, options.altloc);

    const stdout = std.Io.File.stdout();
    var buffer: [8192]u8 = undefined;
    var output = stdout.writer(io, &buffer);
    const writer = &output.interface;
    if (options.format == .tsv) try writer.print("model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_{s}\tpsi_{s}\tomega_{s}\n", .{ unitName(unit), unitName(unit), unitName(unit) });
    if (options.format == .csv) try writer.print("model,chain,residue_number,insertion_code,residue_name,altloc,phi_{s},psi_{s},omega_{s}\n", .{ unitName(unit), unitName(unit), unitName(unit) });
    if (options.format == .json) try writer.writeAll("[\n");

    for (residue_list, 0..) |res, idx| {
        const selected_alt = selected_alts[idx];
        const n = residueAtom(structure, res, "N", selected_alt);
        const ca = residueAtom(structure, res, "CA", selected_alt);
        const c = residueAtom(structure, res, "C", selected_alt);
        // Microheterogeneous sites may use one author residue number for
        // different component IDs in different conformers. Report the
        // component associated with the selected backbone coordinates rather
        // than the first atom-site row seen for that residue number.
        const residue_name = if (ca) |atom| atom.residue else if (n) |atom| atom.residue else if (c) |atom| atom.residue else res.name;
        var phi: ?f64 = null;
        var psi: ?f64 = null;
        var omega: ?f64 = null;
        if (idx > 0 and peptideLinked(structure, residue_list[idx - 1], res, selected_alts[idx - 1], selected_alt))
            phi = torsion4(residueAtom(structure, residue_list[idx - 1], "C", selected_alts[idx - 1]), n, ca, c);
        if (idx + 1 < residue_list.len and peptideLinked(structure, res, residue_list[idx + 1], selected_alt, selected_alts[idx + 1])) {
            const next_n = residueAtom(structure, residue_list[idx + 1], "N", selected_alts[idx + 1]);
            const next_ca = residueAtom(structure, residue_list[idx + 1], "CA", selected_alts[idx + 1]);
            psi = torsion4(n, ca, c, next_n);
            omega = torsion4(ca, c, next_n, next_ca);
        }

        switch (options.format) {
            .tsv, .csv => {
                const separator: u8 = if (options.format == .tsv) '\t' else ',';
                try writer.print("{d}{c}{s}{c}{s}{c}{s}{c}{s}{c}{s}{c}", .{ model, separator, res.chain, separator, res.seq, separator, res.ins, separator, residue_name, separator, selected_alt, separator });
                try writeOptional(writer, phi, unit, "NA");
                try writer.writeByte(separator);
                try writeOptional(writer, psi, unit, "NA");
                try writer.writeByte(separator);
                try writeOptional(writer, omega, unit, "NA");
                try writer.writeByte('\n');
            },
            .json => {
                if (idx != 0) try writer.writeAll(",\n");
                try writer.print("  {{\"model\":{d},\"chain\":", .{model});
                try writeJsonString(writer, res.chain);
                try writer.writeAll(",\"residue_number\":");
                try writeJsonString(writer, res.seq);
                try writer.writeAll(",\"insertion_code\":");
                try writeJsonString(writer, res.ins);
                try writer.writeAll(",\"residue_name\":");
                try writeJsonString(writer, residue_name);
                try writer.writeAll(",\"altloc\":");
                try writeJsonString(writer, selected_alt);
                try writer.writeAll(",\"phi\":");
                try writeOptional(writer, phi, unit, "null");
                try writer.writeAll(",\"psi\":");
                try writeOptional(writer, psi, unit, "null");
                try writer.writeAll(",\"omega\":");
                try writeOptional(writer, omega, unit, "null");
                try writer.writeAll(",\"unit\":");
                try writeJsonString(writer, unitName(unit));
                try writer.writeByte('}');
            },
        }
    }
    if (options.format == .json) try writer.writeAll("\n]\n");
    try writer.flush();
}

fn reportError(err: anyerror) noreturn {
    // A downstream consumer such as `head` may close a pipe early. Treat that
    // normal Unix pipeline condition as success rather than printing noise.
    if (err == error.WriteFailed) std.process.exit(0);
    const exit_code: u8 = switch (err) {
        error.MissingArgument, error.MissingOptionValue, error.InvalidOptionValue, error.UnknownOption, error.InvalidAtomSpec, error.InvalidUnit => 2,
        error.AtomNotFound, error.ModelNotFound, error.NoResidues, error.NoAtoms, error.UnsupportedFormat, error.TruncatedPdbAtom, error.MissingCifColumn, error.InvalidModel, error.InvalidStructureNumber, error.IncompleteCifRow, error.NonContiguousResidue, error.InvalidCompressedInput => 3,
        error.UndefinedGeometry => 4,
        else => 1,
    };
    std.debug.print("zgeom: error: {s}\n", .{@errorName(err)});
    if (exit_code == 2) usage();
    std.process.exit(exit_code);
}

pub fn main(init: std.process.Init) !void {
    const args_z = init.minimal.args.toSlice(init.arena.allocator()) catch |err| reportError(err);
    const args: []const []const u8 = @ptrCast(args_z);
    if (args.len < 2) {
        usage();
        std.process.exit(2);
    }
    if (std.mem.eql(u8, args[1], "-h") or std.mem.eql(u8, args[1], "--help")) {
        usage();
        return;
    }
    if (std.mem.eql(u8, args[1], "-V") or std.mem.eql(u8, args[1], "--version")) {
        const stdout = std.Io.File.stdout();
        try stdout.writeStreamingAll(init.io, "zgeom " ++ build_options.version ++ "\n");
        return;
    }
    if (std.mem.eql(u8, args[1], "backbone"))
        runBackbone(init.io, init.gpa, args[2..]) catch |err| reportError(err)
    else if (std.mem.eql(u8, args[1], "distance") or std.mem.eql(u8, args[1], "angle") or std.mem.eql(u8, args[1], "dihedral"))
        runScalar(init.io, init.gpa, args[1], args[2..]) catch |err| reportError(err)
    else {
        std.debug.print("zgeom: error: unknown command '{s}'\n", .{args[1]});
        usage();
        std.process.exit(2);
    }
}

test "option parsing" {
    const options = try parseOptions(&.{ "x", "--model", "2", "--format", "json", "--unit", "rad" }, 1);
    try std.testing.expectEqual(@as(?u32, 2), options.model);
    try std.testing.expectEqual(Format.json, options.format);
    try std.testing.expectEqual(Unit.radian, options.unit.?);
}
