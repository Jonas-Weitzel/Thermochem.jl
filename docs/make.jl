using Thermochem
using Documenter

DocMeta.setdocmeta!(Thermochem, :DocTestSetup, :(using Thermochem); recursive=true)

makedocs(;
    modules=[Thermochem],
    authors="Jonas Weitzel",
    sitename="Thermochem.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
