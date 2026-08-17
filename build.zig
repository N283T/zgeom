const std = @import("std");

const version = "0.1.0";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("zgeom", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    const options = b.addOptions();
    options.addOption([]const u8, "version", version);
    const exe = b.addExecutable(.{
        .name = "zgeom",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zgeom", .module = mod },
                .{ .name = "build_options", .module = options.createModule() },
            },
        }),
    });
    b.installArtifact(exe);

    const run_step = b.step("run", "Run zgeom");
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| run_cmd.addArgs(args);
    run_step.dependOn(&run_cmd.step);

    const mod_tests = b.addTest(.{ .root_module = mod });
    const exe_tests = b.addTest(.{ .root_module = exe.root_module });
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&b.addRunArtifact(mod_tests).step);
    test_step.dependOn(&b.addRunArtifact(exe_tests).step);

    // Exercise the installed command-line boundary in addition to the unit
    // tests above.  These checks intentionally use only committed fixtures so
    // `zig build test` remains reproducible without network access or a local
    // structure database.
    const pdb_fixture = b.path("tests/fixtures/edge_cases.pdb");
    const cif_fixture = b.path("tests/fixtures/simple.cif");
    const ter_fixture = b.path("tests/fixtures/ter_break.pdb");
    const coherent_altloc_fixture = b.path("tests/fixtures/coherent_altloc.pdb");
    const noncontiguous_fixture = b.path("tests/fixtures/noncontiguous_residue.cif");
    const modified_protein_fixture = b.path("tests/fixtures/modified_protein.pdb");
    const corrupt_gzip_fixture = b.path("tests/fixtures/corrupt.cif.gz");

    const distance_tsv = b.addRunArtifact(exe);
    distance_tsv.addArg("distance");
    distance_tsv.addFileArg(pdb_fixture);
    distance_tsv.addArgs(&.{ "A:10:N", "A:10:CA" });
    distance_tsv.expectStdOutEqual(
        "kind\tmodel\tatom1\tatom2\tvalue\tunit\n" ++
            "distance\t1\tA:10:N\tA:10:CA@B\t1.414214\tangstrom\n",
    );
    test_step.dependOn(&distance_tsv.step);

    const altloc_csv = b.addRunArtifact(exe);
    altloc_csv.addArg("distance");
    altloc_csv.addFileArg(pdb_fixture);
    altloc_csv.addArgs(&.{ "A:10:N", "A:10:CA", "--altloc", "A", "--format", "csv" });
    altloc_csv.expectStdOutEqual(
        "kind,model,atom1,atom2,value,unit\n" ++
            "distance,1,A:10:N,A:10:CA@A,1.000000,angstrom\n",
    );
    test_step.dependOn(&altloc_csv.step);

    const model_json = b.addRunArtifact(exe);
    model_json.addArg("distance");
    model_json.addFileArg(pdb_fixture);
    model_json.addArgs(&.{ "A:10:N", "A:10:CA", "--model", "2", "--format", "json" });
    model_json.expectStdOutEqual(
        "{\"kind\":\"distance\",\"model\":2,\"atoms\":[\"A:10:N\",\"A:10:CA\"],\"value\":1.000000,\"unit\":\"angstrom\"}\n",
    );
    test_step.dependOn(&model_json.step);

    const angle_json = b.addRunArtifact(exe);
    angle_json.addArg("angle");
    angle_json.addFileArg(cif_fixture);
    angle_json.addArgs(&.{ "A:1:N", "A:1:CA", "A:1:C", "--format", "json" });
    angle_json.expectStdOutEqual(
        "{\"kind\":\"angle\",\"model\":1,\"atoms\":[\"A:1:N\",\"A:1:CA\",\"A:1:C\"],\"value\":90.000000,\"unit\":\"degree\"}\n",
    );
    test_step.dependOn(&angle_json.step);

    const dihedral_tsv = b.addRunArtifact(exe);
    dihedral_tsv.addArg("dihedral");
    dihedral_tsv.addFileArg(cif_fixture);
    dihedral_tsv.addArgs(&.{ "A:1:N", "A:1:CA", "A:1:C", "A:2:N" });
    dihedral_tsv.expectStdOutEqual(
        "kind\tmodel\tatom1\tatom2\tatom3\tatom4\tvalue\tunit\n" ++
            "dihedral\t1\tA:1:N\tA:1:CA\tA:1:C\tA:2:N\t180.000000\tdegree\n",
    );
    test_step.dependOn(&dihedral_tsv.step);

    const backbone_tsv = b.addRunArtifact(exe);
    backbone_tsv.addArg("backbone");
    backbone_tsv.addFileArg(cif_fixture);
    backbone_tsv.expectStdOutEqual(
        "model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_degree\tpsi_degree\tomega_degree\n" ++
            "1\tA\t1\t\tALA\t\tNA\t180.000000\t180.000000\n" ++
            "1\tA\t2\t\tGLY\t\t180.000000\tNA\tNA\n",
    );
    test_step.dependOn(&backbone_tsv.step);

    const backbone_json = b.addRunArtifact(exe);
    backbone_json.addArg("backbone");
    backbone_json.addFileArg(cif_fixture);
    backbone_json.addArgs(&.{ "--format", "json" });
    backbone_json.expectStdOutEqual(
        "[\n" ++
            "  {\"model\":1,\"chain\":\"A\",\"residue_number\":\"1\",\"insertion_code\":\"\",\"residue_name\":\"ALA\",\"altloc\":\"\",\"phi\":null,\"psi\":180.000000,\"omega\":180.000000,\"unit\":\"degree\"},\n" ++
            "  {\"model\":1,\"chain\":\"A\",\"residue_number\":\"2\",\"insertion_code\":\"\",\"residue_name\":\"GLY\",\"altloc\":\"\",\"phi\":180.000000,\"psi\":null,\"omega\":null,\"unit\":\"degree\"}\n" ++
            "]\n",
    );
    test_step.dependOn(&backbone_json.step);

    const ter_backbone = b.addRunArtifact(exe);
    ter_backbone.addArg("backbone");
    ter_backbone.addFileArg(ter_fixture);
    ter_backbone.expectStdOutEqual(
        "model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_degree\tpsi_degree\tomega_degree\n" ++
            "1\tA\t1\t\tALA\t\tNA\tNA\tNA\n" ++
            "1\tA\t2\t\tGLY\t\tNA\tNA\tNA\n",
    );
    test_step.dependOn(&ter_backbone.step);

    const coherent_altloc_backbone = b.addRunArtifact(exe);
    coherent_altloc_backbone.addArg("backbone");
    coherent_altloc_backbone.addFileArg(coherent_altloc_fixture);
    coherent_altloc_backbone.expectStdOutEqual(
        "model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_degree\tpsi_degree\tomega_degree\n" ++
            "1\tA\t1\t\tALA\tA\tNA\t180.000000\t180.000000\n" ++
            "1\tA\t2\t\tGLY\t\t180.000000\tNA\tNA\n",
    );
    test_step.dependOn(&coherent_altloc_backbone.step);

    const requested_altloc_backbone = b.addRunArtifact(exe);
    requested_altloc_backbone.addArg("backbone");
    requested_altloc_backbone.addFileArg(coherent_altloc_fixture);
    requested_altloc_backbone.addArgs(&.{ "--altloc", "B" });
    requested_altloc_backbone.expectStdOutEqual(
        "model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_degree\tpsi_degree\tomega_degree\n" ++
            "1\tA\t1\t\tALA\tB\tNA\t-135.000000\t180.000000\n" ++
            "1\tA\t2\t\tGLY\t\t135.000000\tNA\tNA\n",
    );
    test_step.dependOn(&requested_altloc_backbone.step);

    const noncontiguous_backbone = b.addRunArtifact(exe);
    noncontiguous_backbone.addArg("backbone");
    noncontiguous_backbone.addFileArg(noncontiguous_fixture);
    noncontiguous_backbone.expectExitCode(3);
    noncontiguous_backbone.expectStdOutEqual("");
    test_step.dependOn(&noncontiguous_backbone.step);

    const modified_protein_backbone = b.addRunArtifact(exe);
    modified_protein_backbone.addArg("backbone");
    modified_protein_backbone.addFileArg(modified_protein_fixture);
    modified_protein_backbone.expectStdOutEqual(
        "model\tchain\tresidue_number\tinsertion_code\tresidue_name\taltloc\tphi_degree\tpsi_degree\tomega_degree\n" ++
            "1\tA\t1\t\tMSE\t\tNA\t180.000000\t180.000000\n" ++
            "1\tA\t2\t\tGLY\t\t180.000000\tNA\tNA\n",
    );
    test_step.dependOn(&modified_protein_backbone.step);

    const corrupt_gzip = b.addRunArtifact(exe);
    corrupt_gzip.addArg("backbone");
    corrupt_gzip.addFileArg(corrupt_gzip_fixture);
    corrupt_gzip.expectExitCode(3);
    corrupt_gzip.expectStdOutEqual("");
    test_step.dependOn(&corrupt_gzip.step);

    const cli_error = b.addRunArtifact(exe);
    cli_error.addArg("not-a-command");
    cli_error.expectExitCode(2);
    cli_error.expectStdOutEqual("");
    test_step.dependOn(&cli_error.step);

    const missing_atom = b.addRunArtifact(exe);
    missing_atom.addArg("distance");
    missing_atom.addFileArg(cif_fixture);
    missing_atom.addArgs(&.{ "A:999:N", "A:1:CA" });
    missing_atom.expectExitCode(3);
    missing_atom.expectStdOutEqual("");
    test_step.dependOn(&missing_atom.step);

    const undefined_geometry = b.addRunArtifact(exe);
    undefined_geometry.addArg("angle");
    undefined_geometry.addFileArg(cif_fixture);
    undefined_geometry.addArgs(&.{ "A:1:CA", "A:1:CA", "A:1:C" });
    undefined_geometry.expectExitCode(4);
    undefined_geometry.expectStdOutEqual("");
    test_step.dependOn(&undefined_geometry.step);

    // Optional external-oracle regression. PyMOL is intentionally not a
    // project dependency, so this lives behind its own explicit build step.
    const oracle_step = b.step("oracle-test", "Compare 1CRN backbone torsions with installed PyMOL");
    const pymol_oracle = b.addSystemCommand(&.{ "pymol", "-cq", "-r" });
    pymol_oracle.setEnvironmentVariable("ZGEOM_ORACLE_ROOT", b.pathFromRoot("."));
    pymol_oracle.addFileArg(b.path("tests/oracle/compare_1crn_pymol.py"));
    pymol_oracle.step.dependOn(b.getInstallStep());
    oracle_step.dependOn(&pymol_oracle.step);
}
