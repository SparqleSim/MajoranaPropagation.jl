using IJulia

examples_dir = @__DIR__
docs_examples_dir = normpath(joinpath(examples_dir, "..", "docs", "src", "examples"))
template_path = joinpath(examples_dir, "markdown_template.tpl")
kernel_name = "majoranapropagation-docs"
timeout = 600
python = get(ENV, "PYTHON", "python3")

# Use the custom Julia kernel for Julia notebooks to ensure they use the correct project environment.
IJulia.installkernel(
    kernel_name,
    "--project=$(examples_dir)";
    specname = kernel_name,
    displayname = kernel_name,
)

rm(docs_examples_dir; recursive = true, force = true)
mkpath(docs_examples_dir)

# Only these notebooks are executed and included in the documentation.
# Data-heavy or work-in-progress notebooks (e.g. quantum chemistry examples
# requiring FCIDUMP files) are deliberately excluded.
included_notebooks = [
    "0_Majorana-operators.ipynb",
    "1_custom_Majorana_Operators.ipynb",
    "2_fermionic_gates.ipynb",
    "3_trotter_consistent.ipynb",
    "Hubbard_1d.ipynb",
    "Hubbard_2d.ipynb",
    "ad_example.ipynb",
    "imaginary-time.ipynb",
]

notebooks = [joinpath(examples_dir, name) for name in included_notebooks]

sem = Base.Semaphore(5)

@sync for notebook in notebooks
    @async begin
        Base.acquire(sem)
        try
            # Python notebooks should use their default kernel.
            is_python = occursin("\"language\": \"python\"", read(notebook, String))
            kernel_arg = is_python ? `` : `--ExecutePreprocessor.kernel_name=$kernel_name`

            run(`$python -m nbconvert --to markdown \
                --execute \
                $kernel_arg \
                --ExecutePreprocessor.timeout=$timeout \
                --output-dir $docs_examples_dir \
                --template-file $examples_dir/markdown_template.tpl \
                --NbConvertBase.display_data_priority "['image/svg+xml', 'image/png', \
                    'image/jpeg', 'text/markdown', 'text/plain']" \
                $notebook`)
        catch e
            @error "$notebook failed to run with the error: \n $e"
        finally
            Base.release(sem)
        end
    end
end
