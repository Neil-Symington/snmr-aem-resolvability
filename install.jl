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
#   conda writes its launcher's shebang one of two ways, and the choice is made
#   by PATH LENGTH, not by the operating system: if `#!<abs path to conda's own
#   python>` would exceed conda's 127-character limit it writes
#   `#!/usr/bin/env python` instead, and `env` then resolves `python` through
#   PATH. Measured boundary, same machine and Manifest either side: a depot
#   python path of 125 characters gets the absolute shebang and works with NO
#   python on PATH at all; 126 and above gets `env python` and fails.
#   A short depot -- ~/.julia is about 57 characters -- therefore never hits
#   this, which is why GitHub's runners install cleanly and a deep project-local
#   or scratch depot does not.
#   => conda exits 127 before it does anything, the HiQGA build fails, and
#      instantiate takes the whole environment down with it.
#
#   AND A SYSTEM `python` DOES NOT FIX IT. With `env python` resolving to
#   /usr/bin/python3, conda fails differently: ModuleNotFoundError: No module
#   named 'conda'. Only CONDA'S OWN python can run conda. So what step 3 does
#   is not "supply a python" -- it is "put conda's python where env will find
#   it". An earlier version of this comment blamed Debian and Ubuntu having no
#   system python; that was a wrong explanation of a real failure, and an
#   independent review measured the actual cause.
#
# PUTTING CONDA'S BIN ON PATH IS THE FIX. An earlier version of this file
# credited the alias it creates in step 2; that was wrong. On Miniconda 26.9.4
# the installer already ships the alias and step 2 reports "already present",
# yet the build still fails without the PATH change. The alias step is kept as
# a fallback for older Miniconda that may genuinely have lacked it.
# Verified from an empty depot.
# STEP 0: refuse Julia before 1.12. [compat] julia = "1.12" does NOT enforce
# this -- Pkg prints a warning for the active project and carries on, which a
# review proved by installing on 1.11.7 and crashing deep inside
# PrecompileTools with an UndefVarError for Base.StaticData. The pinned
# Manifest resolves PrecompileTools 1.3.4, whose src/invalidations.jl wants
# Base.ReinferUtils or Base.StaticData, and neither exists before 1.12. So the
# check has to be here, where it can actually stop.
if VERSION < v"1.12"
    println("""
    ERROR: this project needs Julia 1.12 or newer; you are on $VERSION.

      Why: the pinned Manifest.toml resolves PrecompileTools 1.3.4, which
      needs Base.ReinferUtils or Base.StaticData. Neither exists before
      1.12, so the install fails inside PrecompileTools rather than here.

      Get 1.12 with juliaup:  juliaup add 1.12 && juliaup default 1.12
      or from https://julialang.org/downloads/
    """)
    exit(1)
end

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
