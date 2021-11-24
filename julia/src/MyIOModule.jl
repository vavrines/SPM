module MyIOModule

using DelimitedFiles

export readtodict
export niceprint
export nicewrite
export replacelineswithstring

function readtodict(filename::String, allowed)
    println("")
    f = open(filename)
    vars = Dict{String,Any}()
    println("Reading config from $filename")
    for line in eachline(f)

        # skip comments
        if length(line) == 0 || line[1] == '#'
            #println("skip comment line")
            continue
        end
        print("\t")
        println(line)

        var, val = split(line, "=")
        stripped = strip(var)
        if stripped in allowed
            vars[stripped] = parse(Float64, val)
        end
    end
    println("")
    return vars
end





function niceprint(A)
    println("")
    #Base.showarray(STDOUT,A,false)
    println(A)
    println("")
end

function nicewrite(filename::String, A)
    writedlm(filename, A, ", ")
end


function replacelineswithstring(
    filename::String,
    newfilename::String,
    keywords::Array{String},
    replaceby::Array{String},
)
    # replaces every line,  that contains the "keyword" with a line that consists of "replaceby"
    # result is written into a new file called "newfilename"


    # file from which to read
    f = open(filename)
    lines = readlines(f)

    # file to which we write
    open(newfilename, "w") do tmpf
        for l in lines
            containssomething = false
            matchid = 0
            for k = 1:size(keywords, 1)
                if occursin(keywords[k], l)
                    matchid = k
                    containssomething = true
                end
            end
            if !(containssomething)
                write(tmpf, l * "\n")
            else
                write(tmpf, replaceby[matchid] * "\n")
            end
        end
    end

end



end
