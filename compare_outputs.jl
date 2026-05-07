using SHA

function main()
    # Compare all VTK files
    kernel_dir = "src/output"
    verif_dir = "verification/output"

    kernel_files = sort(filter(f -> startswith(f, "LinDiffusion") && endswith(f, ".vtk"), readdir(kernel_dir)))
    verif_files = sort(filter(f -> startswith(f, "LinDiffusion") && endswith(f, ".vtk"), readdir(verif_dir)))

    println("Comparing VTK files:")
    println("Kernel: $(length(kernel_files)) files")
    println("Verification: $(length(verif_files)) files")

    # Check how many match exactly
    matches_count = 0
    for (kf, vf) in zip(kernel_files, verif_files)
        kh = bytes2hex(sha256(open(read, joinpath(kernel_dir, kf))))
        vh = bytes2hex(sha256(open(read, joinpath(verif_dir, vf))))
        if kh == vh
            matches_count += 1
        end
    end

    println("\nFile comparison:")
    println("  Matching SHA256: $(matches_count)/$(length(kernel_files))")
    println("  Different: $(length(kernel_files) - matches_count)/$(length(kernel_files))")

    # Sample first file
    kf0 = kernel_files[1]
    vf0 = verif_files[1]
    println("\nFirst file comparison (t=0):")
    k0h = bytes2hex(sha256(open(read, joinpath(kernel_dir, kf0))))
    v0h = bytes2hex(sha256(open(read, joinpath(verif_dir, vf0))))
    println("  Kernel: $k0h")
    println("  Verif:  $v0h")
    println("  Match: $(k0h == v0h)")
end

main()
