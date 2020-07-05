#!/usr/bin/env julia

using ArgParse
using Printf
using BioStructures

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Superimpose a set of protein structures and report a RSMD matrix, in CSV and Mega-compatible formats."
    s.version = "1.1"
    s.add_version = true

    @add_arg_table! s begin
        "--output", "-o"
            help = "output files base name"
            default = "rmsd_matrix"

        "structure"
            help = "PDB, MMCIF or MMTF structural file"
            required = true
            nargs = '*'

    end

    return parse_args(s)
end

function write_csv(filename::String, struct_names::Array{String, 1}, mat::Array{Float64,2})
    open(filename, "w") do output_file
        write(output_file, "structures")
        for struct_name in struct_names
            write(output_file, ",$struct_name")
        end
        write(output_file, "\n")
        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                if j == 1
                    write(output_file, "$(struct_names[i])")
                end
                if i > j
                    write(output_file, ",$(mat[i,j])")
                else
                    write(output_file, ",")
                end
            end
            write(output_file, "\n")
        end
    end
end

function get_format(filename::String)
    format = nothing

    base_name, ext = splitext(basename(filename))

    if lowercase(ext) == ".pdb" || lowercase(ext) == "ent"
        format = PDB

    elseif lowercase(ext) == ".cif"
        format = MMCIF

    elseif lowercase(ext) == ".mmtf"
        format = MMTF

    end

    format

end

function read_structure(filename::String)
    format = get_format(filename)
    structure = nothing

    if isnothing(format)
        return
    end

    try
        structure = read(filename, format)
    catch
        return
    end

    structure

end

function get_basenames(file_names::Array{String,1})
    base_names::Array{String, 1} = []
    for file_name in file_names
        base_name = splitext(basename(file_name))[1]
        push!(base_names, base_name)
    end
    base_names

end

function write_meg(filename::String, struct_names::Array{String,1}, mat::Array{Float64,2})
    open(filename, "w") do output_file
        write(output_file, "#mega\n")
        write(output_file, "!Title: RMSD matrix;\n")
        write(output_file, "!Format DataType=Distance DataFormat=LowerLeft NTaxa=$(length(struct_names));\n")
        write(output_file, "!Description\n")
        write(output_file, "\tRMSD bewteen structures calculated with rmsdmat\n;\n\n")

        for (i, struct_name) in enumerate(struct_names)
            write(output_file, "[$i] #$struct_name\n")
        end

        write(output_file, "\n[     ")
        for i in 1:length(struct_names)
            @printf(output_file, "%9d", i)
        end
        write(output_file, "  ]\n")

        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                if j == 1
                    #write(output_file, "[$i]   ")
                    @printf(output_file, "[%2d]   ", i)
                end
                if i > j
                    @printf(output_file, "%9.3g", mat[i,j])
                else
                    write(output_file, repeat(" ", 9))
                end
            end
            write(output_file, "\n")
        end
    end
end


function main()
    parsed_args = parse_commandline()
    unique_files::Array{String, 1} = unique(parsed_args["structure"])
    n = length(unique_files)
    if n < 2
        write(stderr, "ERROR: At least two distinct structures are required.")
        exit(1)
    end

    rmsd_mat = zeros(n, n) .- 1.0 # -1.0 is a nonsense rmsd value.

    struct1 = nothing
    struct2 = nothing

    # These cicle has to transverse the lower left corner of the matrix, excluding the diagonal.
    for c in 1:(n - 1)
        struct1 = read_structure(unique_files[c])
        if isnothing(struct1)
            write(stderr, "ERROR: $(unique_files[c]) is not a proper structure file. Exiting...\n")
            exit(1)
        end

        for r in (c + 1):n
            struct2 = read_structure(unique_files[r])
            if isnothing(struct2)
                write(stderr, "ERROR: $(unique_files[r]) is not a proper structure file. Exiting...\n")
                exit(1)
            end

            rmsd_mat[r, c] = rmsd(struct1, struct2)
        end
    end

    struct_names = get_basenames(unique_files)

    write_csv(parsed_args["output"]*".csv", struct_names, rmsd_mat)
    write_meg(parsed_args["output"]*".meg", struct_names, rmsd_mat)
end

main()
