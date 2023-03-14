# A Framework for Categories 

The main idea of TensorCategories is to provide abstract methods for categorical computations. TensorCategories.jl provides many generic types and methods which will work fine on custom categories implementing the interface

# The Basics

For basic functionality of `YourCategory <: Category` together with objects
`YourCategoryObject <: CategoryObject` and morphisms `YourMorphism <: CategoryMorphism` you have to minimally build

```julia
struct YourCategory <: Category
  base_ring::Field
end

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

  - $k$-linear 
  - additive
  - abelian
  - monoidal
  - rigid
  - spherical
  - Krull-Schmidt
  - (multi)ring
  - (multi)tensor
  - (multi)fusion

These structures/properties are handled by providing the corresponding methods for any given category. I.e. we require the following methods:

| Structure/Property | Required Structures/Properties | Additional Methods |
|:------------------ | :----------------------------- | :----------------- |
| Additive           |                                | `direct_sum(::YourObject,::YourObject)::YourObject` <br> `zero(::YourCategory)::YourObject` |
| Linear             |                                | `base_ring(::YourCategory)::Field`|
| Abelian            | Additive, Linear               | `kernel(::YourMorphism)::Tuple{YourObject, YourMorphism}` <br> `cokernel(::YourMorphism)::Tuple{YourObject,YourMorphism}`|