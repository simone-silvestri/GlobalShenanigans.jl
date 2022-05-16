using PackageCompiler
create_sysimage(["GLMakie"], sysimage_path="glmakie.so", precompile_execution_file="viz/precompile_plotting.jl")

# julia --project --sysimage glmakie.so