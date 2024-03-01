Intent: to represent Polynomial Functor Moore machines using Catlab.

To do so, we'll create:
  - A schema for Poly, PolyMaps, and MooreMachine
  - Graphviz visualizations for them
  - Operations defined on them, allowing for coproduct, product, tensor, and tri (composition)
  - A host of examples 



Currently, tri, coproduct, product and tensor work on regular polynomials.  The map is currently broken, but will ultimately have things defined on it.

Graphviz visualizations currently work for polys.



We're creating an implementation of Moore machines as representing by polys, in Algebraic Julia.  This is done by creating various schema categories we'll map between: 'FinPoly', 'FinPolyMap', and 'MooreMachine'.  FinPoly consists of directions, positions, and which position maps to which direction.  FinPolyMap consists of two polys, which positions map to which positions, and for each position's poly it maps to, which directions map back to what.  Moore Machine is a graph, with directed edges between nodes.


We've used GraphViz to display these.  Corollas lend themselves naturally to being displayed.  Each position and direction is a node, and if a direction comes from a particular position, there is an edge connecting it.



We've also defined some operations on polys.  Sum, product and tensor already came from Poly.jl, and from these can defined composition ('tri') and \curlyvee.

We plan on defining these on mappings as well.


