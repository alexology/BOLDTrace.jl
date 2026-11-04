# Start the BOLDTrace web application on http://127.0.0.1:8080.
#
#     julia --threads=2 examples/run_app.jl
#
# also works with `include("examples/run_app.jl")` from the VS Code REPL.
# The usage database and the BOLD checkpoint are created in the current
# directory. Stop with Ctrl-C.
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using BOLDTrace

serve(; database = "boldtrace.sqlite", checkpoint = "bold_checkpoint", port = 8080)
