#=----------------------------------------------------------
    Provide structures for diagrams to define Cones and 
    Limits
----------------------------------------------------------=#


abstract type AbstractDiagram <: Functor end

struct KernelDiagram <: AbstractDiagram end