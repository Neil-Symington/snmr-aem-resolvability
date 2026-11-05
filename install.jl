#!/usr/bin/env julia
# install.jl -- set the environment up, working around an upstream packaging bug.
#
#   julia --project=. install.jl
#
# WHY THIS SCRIPT EXISTS. `Pkg.instantiate()` on its own does not work from a
# clean Julia depot, and the error it gives points at this repository rather
# than at the cause. The chain is:
#
#   HiQGA's deps/build.jl runs Conda.add("matplotlib") unconditionally, and
#   then overwrites ENV["PYTHON"] with its own conda python -- so setting
#   PYTHON yourself does nothing.
#   conda's own launcher begins #!/usr/bin/env python, and `env` resolves
#   `python` through PATH.
#   Miniconda ships that `python` inside its OWN bin, but that bin is not on
#   PATH while the HiQGA build runs, and Debian and Ubuntu have no system
#   `python`, only `python3`.
#   => conda exits 127 before it does anything, the HiQGA build fails, and
#      instantiate takes the whole environment down with it.
#
# PUTTING CONDA'S BIN ON PATH IS THE FIX. An earlier version of this file
# credited the alias it creates in step 2; that was wrong. On Miniconda 26.9.4
# the installer already ships the alias and step 2 reports "already present",
# yet the build still fails without the PATH change. The alias step is kept as
# a fallback for older Miniconda that may genuinely have lacked it.
# Verified from an empty depot.
using Pkg

@info "step 1/4  installing packages (the HiQGA build step is expected to fail here)"
try
    Pkg.instantiate()
    @info "instantiate completed with no workaround needed"
catch
    @info "instantiate stopped at a build step, as expected -- continuing"
end

@info "step 2/4  checking conda's python alias (a fallback; the PATH in step 3 is the fix)"
# Conda is a DIRECT dependency of this project purely so this line works: it is
# a transitive dependency of HiQGA, and `using Conda` failed here when it was
# only transitive.
using Conda
bin = joinpath(Conda.ROOTENV, "bin")
isdir(bin) || error("conda was never bootstrapped at $bin; rerun and read the log")
if isfile(joinpath(bin, "python")) || islink(joinpath(bin, "python"))
    @info "  alias already present"
else
    cand = filter(f -> occursin(r"^python3\.\d+$", f), readdir(bin))
    isempty(cand) && error("no python3.x found in $bin")
    symlink(first(cand), joinpath(bin, "python"))
    @info "  created python -> $(first(cand))"
end

@info "step 3/4  building HiQGA with conda's bin on PATH"
ENV["PATH"] = bin * ":" * ENV["PATH"]
Pkg.build("HiQGA")

@info "step 4/4  finishing and loading everything"
Pkg.instantiate()
using HiQGA, SMRPInversion, SNMRForward, JLD2
@info "environment ready. next:  julia --project=. test/runtests.jl"
