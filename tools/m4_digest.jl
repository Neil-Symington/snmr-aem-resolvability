# m4_digest.jl -- identity of a truth file, by file hash AND by content.
#
#   julia --project=. tools/m4_digest.jl <truth.jld2> <win_lo> <win_hi> <pkg_lo> <pkg_hi>
#
# Prints TOML lines for the M4 manifest. Reads only.
#
# Two machines building the same column write files whose bytes differ, so
# hash(read(path)) cannot confirm a column across machines. The content digest
# is SHA-256 over the raw Float64 bytes of the arrays that define the column,
# in a fixed order, so it depends on the values and nothing else.
using JLD2, SHA, Printf

const KEYS = ("z", "σ_bulk", "θ_SMR", "θ_total", "clay", "φ")

f = ARGS[1]
wlo, whi, plo, phi_ = parse.(Float64, ARGS[2:5])
T = load(f)
ctx = SHA2_256_CTX()
for k in KEYS
    v = Vector{Float64}(T[k])
    update!(ctx, Vector{UInt8}(codeunits(k)))
    update!(ctx, reinterpret(UInt8, v))
end
digest = bytes2hex(SHA.digest!(ctx))
z = T["z"]; s = T["σ_bulk"]; dz = z[2] - z[1]
cond(a, b) = sum(s[(z .>= a) .& (z .< b)]) * dz

println("truth_file          = \"", f, "\"")
println("truth_file_hash     = \"", string(hash(read(f)), base = 16), "\"   # julia hash(read(path)); machine-dependent")
println("truth_file_sha256   = \"", bytes2hex(sha256(read(f))), "\"")
println("content_digest      = \"", digest, "\"   # sha256 over ", join(KEYS, ", "))
@printf("window_m            = [%.2f, %.2f]\n", wlo, whi)
@printf("window_conductance  = %.4f   # S, true\n", cond(wlo, whi))
@printf("package_m           = [%.2f, %.2f]\n", plo, phi_)
@printf("package_conductance = %.4f   # S, true\n", cond(plo, phi_))
