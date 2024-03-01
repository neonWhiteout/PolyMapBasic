using Catlab.CategoricalAlgebra
using Catlab.Graphics.Graphviz

using GraphViz
using Poly

using Poly.FinPolynomials

export ⊲, ⋎, FinPolyMap, GraphPoly


"""
Schema for map between two fin polys.
Made up of two fin polys, a morphism from the source position
to the target position, and a φ♯ object, which has a morphism to
srcpos, tgtpos, and srcdir, representing a map back which takes into account
the original position.
"""
@present FinPolyMapSchema(FreeSchema) begin
    (SrcPos, SrcDir, TgtPos, TgtDir, φ♯)::Ob
    srcpoly::Hom(SrcDir, SrcPos)
    tgtpoly::Hom(TgtDir, TgtPos)
    φ::Hom(SrcPos, TgtPos)
    φ♯sp::Hom(φ♯, SrcPos)
    φ♯td::Hom(φ♯, TgtDir)
    φ♯sd::Hom(φ♯, SrcDir)    
end


@abstract_acset_type AbstractFinPolyMap
@acset_type FinPolyMap(FinPolyMapSchema, index=[:srcpoly, :tgtpoly, :φ, :φ♯]) <: AbstractFinPolyMap

"""
    FinPolyMap(srcdirs::Vector{Int}, tgtdirs::Vector{Int}, srcmap::Vector{Int}, tgtmap::Vector{Pair{Int,Vector{Int}}}))

Constructor for subtypes of AbstractFinPolyMap.  First two arguments take the
format of FinPoly constructors, with an ordered list of exponents.
[2,3,0,0] = y² + y³ + 1 + 1.  srcmap is a list which says what each index of the
list maps to.  That is, first index of source maps to the index of the int here.
tgtmap takes the compressed form of pairs of ints and vectors of ints.
Index corresponds to srcpos, the int in the pair tgtpos, each index in the inner
list is the relative direction on that tgtpos, and each value in the list is the
relative srcdir.

# Examples
```julia-repl
julia> tgtmap = [1 => [2,3], 1=> [1,1], 4 => Vector{Int}()]^C

julia> srca = [3,1,0];

julia> tgta = [2,1,0,0];

julia> srcmap = [1,1,4];

julia> tgtmap = [1 => [2,3], 1=> [1,1], 4 => Vector{Int}()];

julia> ok = FinPolyMap(srca, tgta, srcmap, tgtmap)
FinPolyMap {SrcPos:3, SrcDir:4, TgtPos:4, TgtDir:3, φ♯:4}
┌────────┬───┐
│ SrcPos │ φ │
├────────┼───┤
│      1 │ 1 │
│      2 │ 1 │
│      3 │ 4 │
└────────┴───┘
┌────────┬─────────┐
│ SrcDir │ srcpoly │
├────────┼─────────┤
│      1 │       1 │
│      2 │       1 │
│      3 │       1 │
│      4 │       2 │
└────────┴─────────┘
┌────────┬─────────┐
│ TgtDir │ tgtpoly │
├────────┼─────────┤
│      1 │       1 │
│      2 │       1 │
│      3 │       2 │
└────────┴─────────┘
┌────┬──────┬──────┬──────┐
│ φ♯ │ φ♯sp │ φ♯td │ φ♯sd │
├────┼──────┼──────┼──────┤
│  1 │    1 │    1 │    2 │
│  2 │    1 │    2 │    3 │
│  3 │    2 │    1 │    4 │
│  4 │    2 │    2 │    4 │
└────┴──────┴──────┴──────┘

```

"""
function (::Type{P})(srcdirs::Vector{Int}, tgtdirs::Vector{Int}, srcmap::Vector{Int}, tgtmap::Vector{Pair{Int,Vector{Int}}}) where P <: AbstractFinPolyMap
    srcpos = length(srcdirs)
    tgtpos = length(tgtdirs)

    tgtdir = length(tgtdirs)
    srcdir = length(srcdirs)
    
    p = P()
    add_parts!(p, :TgtPos, tgtpos)

    add_parts!(p, :SrcPos, srcpos; φ=srcmap) 

    for (pos, dir) in enumerate(srcdirs)
        if dir == 0
            continue
        end
        add_parts!(p, :SrcDir, dir; srcpoly=fill(pos,dir))
        end
      s = 0
    for (pos, dir) in enumerate(tgtdirs)
        add_parts!(p, :TgtDir, dir; tgtpoly=pos)
    end

    for srcpos in 1:length(tgtmap)
        srcdir = incident(p, srcpos, :srcpoly)
        tgtpos = first(tgtmap[srcpos])
        tgtdir = incident(p, tgtpos, :tgtpoly)
        direction_maps = last(tgtmap[srcpos])
        for tgtdirindex in 1:length(direction_maps)
            globaltgtdir = tgtdir[tgtdirindex]
            globalsrcdir = srcdir[direction_maps[tgtdirindex]]
            add_part!(p, :φ♯ ; φ♯sp = srcpos, φ♯td = globaltgtdir, φ♯sd = globalsrcdir)
        end
    end
    
  p
end





srcpos(m::AbstractFinPolyMap) = nparts(m, :SrcPos)
tgtpos(m::AbstractFinPolyMap) = nparts(m, :TgtPos)

"""
Given a poly map, create its srcpoly
"""
create_srcpoly(m::AbstractFinPolyMap) = begin
    exponents = [length(incident(m, i, :srcpoly)) for i in 1:(srcpos(m))]
    FinPoly(exponents)
end

"""
Given a poly map, create its tgtpoly
"""
create_tgtpoly(m::AbstractFinPolyMap) = begin
    exponents = [length(incident(m, i, :tgtpoly)) for i in 1:(tgtpos(m))]
    FinPoly(exponents)
end


subscript_replace = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄', '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉');
superscript_replace = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹');


function polyStatements(y::AbstractFinPoly ;charval='y', colours = ["teal", "yellow", "orange"])
    NNodes = [Node("n$charval$n", Attributes(:label=>"$( charval * replace(string(length(incident(y, n, :pos))), superscript_replace...) * "   " *replace(string(n), subscript_replace...))", :shape=>"square", :color => colours[1])) for n in positions(y)]
    Morphs = []
    
    for n in positions(y)
        for p in incident(y, n, :pos)
            push!(NNodes, Node("p$charval$p", Attributes(:label=>"$(string(p) * replace(string(n), subscript_replace...))",:shape=>"square", :color => colours[2])))
            push!(Morphs, Edge(["n$charval$n", "p$charval$p"], Attributes(:color => colours[3])))
        end 
    end
    stmts = vcat(NNodes, Morphs)
    return stmts
end

"""
Graphs a FinPoly.  Optional keyword argument charval to change symbol label.
"""
function GraphPoly(y::AbstractFinPoly; charval='y')
    stmts = polyStatements(y;charval=charval);
    graph_attrs = Attributes(:rankdir=>"BT")
    edge_attrs  = Attributes(:splines=>"splines")
    return Graphviz.Digraph("G", stmts...; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
end

"""
Composition operator.  Relies on + and * defined in Poly.jl.
"""
⊲(A::AbstractFinPoly, B::AbstractFinPoly) = begin
    sum = FinPoly(Vector{Int}())
    for p in positions(A)
        prod = FinPoly([0])
        for j in incident(A, p, :pos)
            prod = prod * B
        end
        sum = sum + prod
    end
    return sum
end

"""
Vee operator.  Relies on + and ⊗ in Poly.jl.
"""
⋎(A::AbstractFinPoly, B::AbstractFinPoly) = A + A ⊗ B + B


