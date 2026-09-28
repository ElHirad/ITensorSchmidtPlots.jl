using TOML
length(ARGS) == 3 || error("Usage: julia security/audit_dependencies.jl Manifest.toml /path/to/SecurityAdvisories.jl report.toml")
manifest = TOML.parsefile(ARGS[1])
root = ARGS[2]
versions = Dict(k => get(only(v),"version",nothing) for (k,v) in manifest["deps"])
operators = Dict("<" => (<), "<=" => (<=), "=" => (==), ">=" => (>=), ">" => (>))
function affected(version, range)
    strip(range) == "*" && return true
    clauses = map(split(range,',')) do part
        m = match(r"^\s*(<=|>=|<|>|=)\s*(\S+)\s*$",part)
        isnothing(m) && error("Unsupported advisory range: $range")
        (operators[m[1]], VersionNumber(m[2]))
    end
    return all(f(version,bound) for (f,bound) in clauses)
end
@assert affected(v"1.0.0+0", ">= 1.0.0+0, < 1.0.0+1")
@assert !affected(v"1.0.0+1", ">= 1.0.0+0, < 1.0.0+1")
@assert affected(v"1.0.0", "<= 1.0.0")
@assert !affected(v"1.0.0", "< 1.0.0")
findings = Dict{String,Any}[]
counts = Dict("advisories"=>0,"candidate_ranges"=>0,"matched_ranges"=>0,"versioned_packages"=>count(!isnothing,values(versions)))
for (dir,_,files) in walkdir(joinpath(root,"advisories","published")), file in files
    endswith(file,".md") || continue
    path=joinpath(dir,file)
    source=read(path,String)
    front=match(r"(?ms)^`{3,4}toml\s*\n(.*?)\n`{3,4}",source)
    isnothing(front) && error("Missing advisory metadata: $file")
    adv=TOML.parse(front[1])
    counts["advisories"]+=1
    haskey(adv,"withdrawn") && continue
    for item in get(adv,"affected",[])
        pkg=item["pkg"]
        installed=get(versions,pkg,nothing)
        isnothing(installed) && continue
        for range in item["ranges"]
            counts["candidate_ranges"]+=1
            if affected(VersionNumber(installed),range)
                counts["matched_ranges"]+=1
                title=match(r"(?m)^# (.+)$",source)
                push!(findings,Dict("id"=>adv["id"],"package"=>pkg,"installed"=>installed,
                    "range"=>range,"aliases"=>get(adv,"aliases",String[]),"upstream"=>get(adv,"upstream",String[]),
                    "summary"=>isnothing(title) ? "" : title[1],
                    "source"=>"https://github.com/JuliaLang/SecurityAdvisories.jl/blob/main/" * replace(relpath(path,root), '\\'=>'/')))
            end
        end
    end
end
open(ARGS[3],"w") do io
    TOML.print(io,Dict("counts"=>counts,"findings"=>findings); sorted=true)
end
println(counts)
println("These are version-range matches, not proof of reachability or exploitability; see docs/SECURITY_REVIEW.md.")
for f in findings
    println(f["package"]," ",f["installed"]," | ",f["id"]," | ",f["range"]," | ",f["summary"])
end
