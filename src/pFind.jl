module pFind

import CSV
import DataFrames
import MesMS
import RelocatableFolders: @path

const DIR_DATA = @path joinpath(@__DIR__, "../data/pFind")

read_mgf(io::IO) = begin
    M = MesMS.MS2[]
    id = 0
    mz = 0.0
    z = 0
    peaks = MesMS.Peak[]
    while !eof(io)
        line = readline(io)
        if length(line) == 0
            continue
        elseif line == "BEGIN IONS"
            id = 0
            mz = 0.0
            z = 0
            peaks = MesMS.Peak[]
        elseif line == "END IONS"
            push!(M, MesMS.MS2(; id, ions=[MesMS.Ion(mz, z)], peaks))
        elseif startswith(line, "TITLE=")
            id = parse(Int, split(line[7:end], '.')[end-4])
        elseif startswith(line, "PEPMASS=")
            mz = parse(Float64, line[9:end])
        elseif startswith(line, "CHARGE=")
            z = line[end] == '+' ? parse(Int, line[8:end-1]) : -parse(Int, line[8:end-1])
        elseif !occursin('=', line)
            m, i = split(line)
            push!(peaks, MesMS.Peak(parse(Float64, m), parse(Float64, i)))
        end
    end
    return M
end

read_mgf(fname::AbstractString) = open(read_mgf, fname)

read_element(path=joinpath(DIR_DATA, "element.ini")) = begin
    @info "Element loading from " * path
    lines = map(strip, readlines(path)[begin+1:end])
    filter!(l -> startswith(l, "E"), lines)
    d = map(lines) do line
        attrs = split(strip(split(line, '='; limit=2)[2]), '|')[begin:end-1]
        name = Symbol(attrs[1])
        mass = map(s -> parse(Float64, s), split(attrs[2], ',')[begin:end-1])
        abu = map(s -> parse(Float64, s), split(attrs[3], ',')[begin:end-1])
        return name => mass[argmax(abu)]
    end
    return Dict(d)
end

read_amino_acid(path=joinpath(DIR_DATA, "aa.ini")) = begin
    @info "AA loading from " * path
    lines = map(strip, readlines(path)[begin+1:end])
    filter!(l -> startswith(l, "R"), lines)
    d = map(lines) do line
        attrs = split(strip(split(line, '='; limit=2)[2]), '|')[begin:end-1]
        return Symbol(attrs[1]) => parse(MesMS.Formula, attrs[2])
    end
    return Dict(d)
end

read_mod(path=joinpath(DIR_DATA, "modification.ini")) = begin
    @info "Mod. loading from " * path
    lines = map(strip, readlines(path)[begin+1:end])
    filter!(l -> !startswith(l, "name") && !isempty(l), lines)
    d = map(lines) do line
        name, attrs = map(strip, split(line, '='; limit=2))
        attrs = split(attrs)
        site = attrs[1]
        pos = attrs[2]
        mass = parse(Float64, attrs[3])
        comp = parse(MesMS.Formula, attrs[end])
        return Symbol(name) => (; name, site, pos, mass, comp)
    end
    return Dict(d)
end

parse_title(title) = begin
    file, scan, _, _, idx, _  = rsplit(strip(title), '.'; limit=6)
    return file, parse(Int, scan), parse(Int, idx)
end

parse_mod(pep, mod) = begin
    mods = map(m -> match(r"^(\d+),(.+)$", String(m)).captures, split(mod === missing ? "" : mod, ';'; keepempty=false))
    return sort!(map(m -> (Symbol(m[2]), max(min(length(pep), parse(Int, m[1])), 1)), mods))
end

read_psm(path; silencewarnings=false) = begin
    @info "PSM loading from " * path
    df = DataFrames.DataFrame(CSV.File(path; delim='\t', missingstring=nothing, silencewarnings))
    DataFrames.rename!(df, Dict(
        :File_Name => :title,
        :Charge => :z,
        Symbol("Exp.MH+") => :mh,
        Symbol("Calc.MH+") => :mh_calc,
        Symbol("Q-value") => :q_value,
        :Raw_Score => :score_raw,
        :Final_Score => :score,
        :Proteins => :prot,
        Symbol("Target/Decoy") => :td,
    ))
    df.id = Vector{Int}(1:DataFrames.nrow(df))
    DataFrames.transform!(df, :title => DataFrames.ByRow(parse_title) => [:file, :scan, :idx_pre])
    DataFrames.transform!(df, :Sequence => DataFrames.ByRow(MesMS.unify_aa_seq) => :pep)
    DataFrames.transform!(df, [:pep, :Modification] => DataFrames.ByRow(parse_mod) => :mod)
    DataFrames.transform!(df, [:mh_calc, :mh] => DataFrames.ByRow(MesMS.error_ppm) => :error)
    DataFrames.transform!(df, [:mh, :z] => DataFrames.ByRow(MesMS.mh_to_mz) => :mz)
    DataFrames.transform!(df, [:mh_calc, :z] => DataFrames.ByRow(MesMS.mh_to_mz) => :mz_calc)
    df.td = Symbol.(df.td)
    return df
end

modstr(mods; null="") = isempty(mods) ? null : join(map(m -> "$(m[1])@$(m[2])", mods), ",")

pepstr(seq, mods; delim="") = begin
    if isempty(mods)
        return seq
    else
        return "$(seq)$(delim)($(modstr(mods)))"
    end
end

find_prot(seq, fasta) = begin
    prot = [k for (k, v) in fasta if occursin(seq, v)]
    return (; prot=join(prot, '/'), td=!any(p -> startswith(p, "REV_"), prot))
end

find_prot!(seq, fasta, mem) = begin
    if !haskey(mem, seq)
        mem[seq] = find_prot(seq, fasta)
    end
    return mem[seq]
end

parse_top_n_row(line; itol=true) = begin
    items = split(line, '\t')
    pep = items[2]
    n_mod = parse(Int, items[7])
    mods = map(1:n_mod) do i
        mod = Symbol(rsplit(items[7 + 2 * i], '#'; limit=2)[1])
        site = min(length(pep), max(1, parse(Int, items[7 + 2 * i - 1])))
        return (mod, site)
    end
    return (; pep=MesMS.unify_aa_seq(pep; itol), mod=mods)
end

read_top_n(path, fasta=nothing; itol=true) = begin
    lines = open(path) do io
        return collect(eachline(io))
    end
    D = typeof((; file="", scan=0, idx_pre=0, mz=0.0, z=0, iden=[]))[]
    i = 1
    mem = Dict()
    while i ≤ length(lines)
        if !startswith(lines[i], "S\t")
            throw("unexcepted line: $(lines[i])")
        end
        items = split(lines[i], '\t')
        mz = parse(Float64, items[2])
        z = parse(Int, items[3])
        i += 1
        title = strip(lines[i])
        file, scan, idx_pre = pFind.parse_title(title)
        i += 1
        iden = []
        while i ≤ length(lines) && !startswith(lines[i], "S\t")
            push!(iden, lines[i])
            i += 1
        end
        iden = parse_top_n_row.(iden; itol)
        if !isnothing(fasta)
            iden = [(; i..., find_prot!(i.pep, fasta, mem)...) for i in iden]
        end
        push!(D, (; file, scan, idx_pre, mz, z, iden))
    end
    return D
end

end
