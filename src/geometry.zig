const std = @import("std");

pub const Vec3 = struct { x: f64, y: f64, z: f64 };

fn sub(a: Vec3, b: Vec3) Vec3 {
    return .{ .x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z };
}

fn dot(a: Vec3, b: Vec3) f64 {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

fn cross(a: Vec3, b: Vec3) Vec3 {
    return .{
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x,
    };
}

fn norm(a: Vec3) f64 {
    return @sqrt(dot(a, a));
}

pub fn distance(a: Vec3, b: Vec3) f64 {
    return norm(sub(a, b));
}

pub fn angle(a: Vec3, b: Vec3, c: Vec3) ?f64 {
    const ba = sub(a, b);
    const bc = sub(c, b);
    const denominator = norm(ba) * norm(bc);
    if (denominator <= std.math.floatEps(f64)) return null;
    return std.math.acos(std.math.clamp(dot(ba, bc) / denominator, -1.0, 1.0));
}

/// Signed IUPAC torsion in [-pi, pi]. Positive means the clockwise rotation of
/// the far bond when viewed from atom 2 toward atom 3.
pub fn dihedral(a: Vec3, b: Vec3, c: Vec3, d: Vec3) ?f64 {
    const b1 = sub(b, a);
    const b2 = sub(c, b);
    const b3 = sub(d, c);
    const n1 = cross(b1, b2);
    const n2 = cross(b2, b3);
    const b1_len = norm(b1);
    const b2_len = norm(b2);
    const b3_len = norm(b3);
    if (b1_len <= std.math.floatEps(f64) or b2_len <= std.math.floatEps(f64) or b3_len <= std.math.floatEps(f64)) return null;
    // Test collinearity relative to the bond-vector scale. An absolute test on
    // cross-product magnitude would classify identical geometry differently
    // after a coordinate-unit change.
    const collinear_tolerance = 1e-12;
    if (norm(n1) / (b1_len * b2_len) <= collinear_tolerance or
        norm(n2) / (b2_len * b3_len) <= collinear_tolerance) return null;
    const b2n = Vec3{ .x = b2.x / b2_len, .y = b2.y / b2_len, .z = b2.z / b2_len };
    const m = cross(b2n, n1);
    return std.math.atan2(dot(m, n2), dot(n1, n2));
}

pub fn radiansToDegrees(value: f64) f64 {
    return value * 180.0 / std.math.pi;
}

test "geometry definitions" {
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), distance(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 3, .y = 4, .z = 0 }), 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, std.math.pi / 2.0), angle(.{ .x = 1, .y = 0, .z = 0 }, .{ .x = 0, .y = 0, .z = 0 }, .{ .x = 0, .y = 1, .z = 0 }).?, 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, std.math.pi / 2.0), dihedral(.{ .x = 0, .y = 1, .z = 0 }, .{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1, .y = 0, .z = 0 }, .{ .x = 1, .y = 0, .z = 1 }).?, 1e-12);
    try std.testing.expect(angle(.{ .x = 0, .y = 0, .z = 0 }, .{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1, .y = 0, .z = 0 }) == null);
}

test "dihedral degeneracy classification is scale independent" {
    const small = dihedral(.{ .x = 0, .y = 1e-6, .z = 0 }, .{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1e-6, .y = 0, .z = 0 }, .{ .x = 1e-6, .y = 0, .z = 1e-6 });
    const large = dihedral(.{ .x = 0, .y = 1e6, .z = 0 }, .{ .x = 0, .y = 0, .z = 0 }, .{ .x = 1e6, .y = 0, .z = 0 }, .{ .x = 1e6, .y = 0, .z = 1e6 });
    try std.testing.expect(small != null and large != null);
    try std.testing.expectApproxEqAbs(small.?, large.?, 1e-12);
}
