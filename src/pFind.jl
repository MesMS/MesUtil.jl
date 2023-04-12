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
    raw, scan, _, _, idx, _  = rsplit(strip(title), '.'; limit=6)
    return raw, parse(Int, scan), parse(Int, idx)
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
        :Proteins => :proteins,
        Symbol("Target/Decoy") => :td,
    ))
    df.id = Vector{Int}(1:DataFrames.nrow(df))
    DataFrames.transform!(df, :title => DataFrames.ByRow(parse_title) => [:raw, :scan, :idx_pre])
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

end
