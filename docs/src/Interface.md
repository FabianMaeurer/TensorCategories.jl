# A Framework for Categories 

The main incentive of TensorCategories.jl is to provide abstract methods for categorical computations. TensorCategories.jl provides many generic types and methods which will work fine on custom categories implementing the interface.

# The Basics

We provide the following abstract types to identify categories, objects and morphisms:

```julia
abstract type Category end
abstract type CategoryObject end
abstract type CategoryMorphism end
```

For basic functionality of `YourCategory <: Category` together with objects
`YourCategoryObject <: CategoryObject` and morphisms `YourMorphism <: CategoryMorphism` you have to minimally build

```julia
struct YourCategory <: Category end

struct YourCategoryObject <: CategoryObject
  parent::YourCategory
end

struct YourMorphism <: CategoryMorphism
  domain::YourCategoryObject
  codomain::YourCategoryObject
end
```

For objects you need to provide the following methods:

```
id(X::YourCategoryObject) ::YourMorphism
```

For morphisms you need

```
compose(f::YourMorphism, g::YourMorphism) ::YourMorphism
```

# Categories with Additional Structure

At the current state TensorCategories.jl knows about the following kinds of categories:

  - ``k``-linear 
  - additive
  - abelian
  - monoidal
  - rigid monoidal
  - braided monoidal
  - spherical
  - Krull-Schmidt
  - (multi)ring
  - (multi)tensor
  - (multi)fusion

For the precise definitions take a look at [EGNO](@ref), which we use as a general reference in most concerns. The above labels are generically tested by functions `is_label(::Category)::Bool`. In general it is only necessary to implement the necessary methods for a specific type. In some special cases it might be useful to explicitly falsify some labels. For example it is not feasible to check whether a preabelian category is in fact abelian. Thus TensorCategories will always see a preabelian category as abelian by omitting to check for (co)normality of epi/monomorphisms. 

| Structure/Property | Required Structures/Properties | Defining Methods   | Useful Methods  | 
|:------------------ | :----------------------------- | :----------------- | :-------------- |
| Additive           |                                | `direct_sum(::YourObject,::YourObject)::Tuple{YourObject, Vector{YourMorphism}, Vector{YourMorphism}`  `direct_sum(::YourMorphism,::YourMoprhism)::YourMorphism`  `zero(::YourCategory)::YourObject`  `+(::YourMorphism, ::YourMorphism)`| |
| Linear             |                                | `base_ring(::YourCategory)::Field`  `*(::FieldElem,::YourMorphism`  `+(::YourMorphism, ::YourMorphism)` | |
| Abelian            | Additive                       | `kernel(::YourMorphism)::Tuple{YourObject, YourMorphism}`  `cokernel(::YourMorphism)::Tuple{YourObject,YourMorphism}`| `indecomposables(::YourCategory)::Vector{YourObject}`  `simples(::YourCategory)::Vector{YourObject}`| 
| Monoidal           |                                | `tensor_product(::YourObject,::YourObject)::YourObject`  `tensor_product(::YourMorphism,::YourMorphism)::YourMorphism`  `one(::YourCategory)::YourObject` | | 
| Rigid              | Monoidal                       | `left_dual(X::YourObject)::YourObject`  `right_dual(X::YourObject)::YourObject`  or, if the dual is two sided  `dual(::YourObject)::YourObject`  `ev(::YourObject)::YourMorphism`  `coev(::YourObject)::YourMorphism` | |
| Braided            | Monoidal                       | `braiding(X::YourObject, Y::YourObject)::YourMorphism`          | |
| Spherical          | Monoidal, Rigid                | `spherical(X::YourObject)::YourMorphism` |
| Krull-Schmidt      | Additive, Linear               | **Coming Soon** ||
| (Multi)Ring        | Abelian, Monoidal, linear      | | |
| (Multi)Tensor      | Multiring, Rigid               | | |
| (Multi)Fusion      | Multitensor, Semisimple        | | |

!!! warn "Beware"
    If any of the above functions are implemented only partially it might be useful to manually specify the according label to be false. Otherwise there might be unwanted behavior. 

!!! note "Direct sums"
    Direct sums always return the object and the corresponding inclusions and projections.