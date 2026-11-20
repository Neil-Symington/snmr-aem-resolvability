# runout.jl -- every stage writes a RUN FOLDER, and no stage overwrites a result.
#
# Ruled by Neil, 21 Nov 2026, after the third artefact lost to an overwrite:
# Stage 2 run A (18 Nov), the control nuisance files (18 Nov), and the 2P
# production gate (21 Nov, overwritten by a three-node smoke that shared its
# config and therefore its output path). The case was fixed twice; this fixes
# the class.
#
# The shape, which Stage 2 already had:
#
#   <out.dir>/runs/<stage>_<name>_<YYYYmmdd-HHMMSS>/   the run, immutable
#   <out.dir>/<stage>_<name>.jld2                      symlink -> newest run
#   <out.dir>/<stage>_<name>_LATEST                    text pointer, one line
#
# The flat path stays, as a symlink, so every reader that already knows it --
# the report pages, the acceptance test, the figure scripts -- keeps working.
# The POINTER moves; the RESULT never does. That is the distinction the ruling
# draws, and the reason there is no --force anywhere in this file.
module RunOut
using Dates
export run_dir, save_result, publish_run

"A fresh, empty run folder under `base`. Never an existing one."
function run_dir(base::AbstractString, stage::AbstractString, name::AbstractString)
    stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
    d = joinpath(base, "runs", "$(stage)_$(name)_$stamp")
    k = 0
    while ispath(d)                      # two launches inside one second
        k += 1
        d = joinpath(base, "runs", "$(stage)_$(name)_$(stamp)_$k")
    end
    mkpath(d)
    d
end

"""
    save_result(f, path)

Call `f(path)` only if nothing is there. A result file is immutable: there is
no force, and a stage that wants to write again writes to a new run folder.
"""
function save_result(f::Function, path::AbstractString)
    ispath(path) && error("""
        refusing to overwrite an existing result:
          $path
        Results are immutable (ruling of 21 Nov 2026). Write into a new run
        folder from RunOut.run_dir instead. There is no --force.""")
    mkpath(dirname(path))
    f(path)
    path
end

"""
    publish_run(base, dir) -> Vector{String}

Point `base` at the run in `dir`: one symlink per file produced, replacing any
earlier symlink, plus a one-line `_LATEST` pointer naming the run folder.
Replacing a symlink is not overwriting a result -- the results stay in their
run folders, and every earlier run remains readable at its own path.
"""
function publish_run(base::AbstractString, dir::AbstractString)
    made = String[]
    for f in sort(readdir(dir))
        src = joinpath(dir, f)
        isfile(src) || continue
        link = joinpath(base, f)
        if islink(link) || isfile(link)
            islink(link) || @warn "replacing a plain file with a pointer; the \
                                   result itself is kept in its run folder" path=link
            rm(link)
        end
        symlink(relpath(src, base), link)
        push!(made, link)
    end
    open(joinpath(base, "LATEST"), "w") do io
        println(io, basename(dir))
    end
    made
end
end
