# Generates HTML documentation from the contents of
# the docs folder. To generate, we must first setup
# a symlink from the repo README.md to src/index.md,
# in order re-use the README in Documenter.jl
# From the docs/ directory (containing this file):
#     cd src
#     ln -s ../../README.md index.md
#     cd ../
#
# This need only be done once per-machine. Then,
# generating/updating the doc is triggered via
#     julia --project make.jl
#
# If triggered within a Github Action, the generated
# HTML files will then be committed to the 'gh-pages'
# branch, which Github Pages can be configured to
# display at SparqleSim.github.io/MajoranaPropagation.jl/
#
# Note documentation generated from non-main branches
# will be uploaded to subdomain /dev/, even when not
# from the 'dev' branch, and doc generated from pull
# requests will be uploaded to /previews/PR#.

import Pkg

docs_dir = @__DIR__
package_dir = normpath(joinpath(docs_dir, ".."))
examples_dir = joinpath(docs_dir, "..", "examples")

function build_examples()
    Pkg.activate(examples_dir)
    Pkg.develop(Pkg.PackageSpec(path = package_dir))
    Pkg.instantiate()

    include(joinpath(examples_dir, "make.jl"))
end

build_examples()


Pkg.activate(docs_dir)
using Documenter, MajoranaPropagation


# Generate doc HTML files, saved to build/
makedocs(
    # Add favicon.ico
    format=Documenter.HTML(
        assets=[
            "assets/favicon.ico",
            "assets/custom.css",
        ],
        size_threshold=nothing,
        size_threshold_warn=200 * 1024,
    ), sitename="MajoranaPropagation.jl",

    # determines site layout
    pages=[

        # index.md does not exist; it is a symlink
        # to the repo's README.md file, created as
        # per the comments above, to avoid duplicating
        # the README.md contents into Documenter.jl
        # pages. We manually override its name in the
        # left navbar to be "Home"
        "Home" => "index.md",

        # these files are automatically created by
        # `build_examples` function using `nbconvert`.
        # the filenames are same as the corresponding
        # .ipynb files. The notebooks have no top-level
        # markdown heading, so navbar titles are set
        # explicitly here.
        "Examples" => [
            "Majorana Operators" => "examples/0_Majorana-operators.md",
            "Custom Majorana Operators" => "examples/1_custom_Majorana_Operators.md",
            "Fermionic Gates" => "examples/2_fermionic_gates.md",
            "Trotterized Time Evolution" => "examples/3_trotter_consistent.md",
            "1D Hubbard Model" => "examples/Hubbard_1d.md",
            "2D Hubbard Model" => "examples/Hubbard_2d.md",
            "Automatic Differentiation" => "examples/ad_example.md",
            "Imaginary-Time Evolution" => "examples/imaginary-time.md",
        ],

        # these 'lower-level' files also exist, and will
        # be grouped under an 'API' section in the navbar
        "API" => [
            "api/MajoranaDataTypes.md",
            "api/MajoranaAlgebra.md",
            "api/InitialStates.md",
            "api/Gates.md",
            "api/Propagation.md",
            "api/Truncations.md",
            "api/Circuits.md",
            "api/FrequencyTracker.md",
            "api/QuantumChemistry.md",
        ]
    ]
)


# When run from a Github Action, commit those files to the 'gh-pages' branch,
# depending upon the triggering branch or whether it is a release/pull-request.
deploydocs(
    repo="github.com/SparqleSim/MajoranaPropagation.jl.git",

    # Enable generation of doc from PRs, under a /previews/PR## sub-domain.
    # Beware that this requires the Github Action was explicitly triggered by
    # a 'pull_request' event (not a 'push')
    push_preview=true,

    # Specify that changes to our 'dev' branch (rather than default 'main')
    # should update the doc visible at the /dev/ sub-domain. Note this means
    # pushes to the main branch never generate doc; only new releases will
    # (see below)
    devbranch="dev",
    devurl="dev",

    # Control which Github releases (here, all) trigger re-generation of the
    # main documentation, and their URLs. Below, we specify that:
    # - subdomain /stable/ presents the very latest release doc
    # - all versions (including patches) have their own hosted doc under /vX.Y.Z/
    # - changes to the dev branch should update the /dev/ sub-domain
    versions=["stable" => "v^", "v#.#.#", "dev" => "dev"]
)


# Once 'gh-pages' branch is updated, and Github Pages has been configured to
# publish files from that branch, the documentation is visible at either:
# - SparqleSim.github.io/MajoranaPropagation.jl/
# - SparqleSim.github.io/MajoranaPropagation.jl/dev/
# - SparqleSim.github.io/MajoranaPropagation.jl/previews/PR#
# where # above is replaced with the pull request number.
#
# These "doc clones" are deleted whenever a commit is pushed to the main
# branch (signifying a version release), so that development history does
# not bloat the repo. Deletion is performed by the 'tidy-doc' CI job.
