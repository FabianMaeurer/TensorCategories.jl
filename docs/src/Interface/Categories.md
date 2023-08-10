```@setup Ex
using TensorCategories, Oscar
```

# Basics

The Interface is naturally based on three abstract types which 
have to be extended:

```julia 
abstract type Category end
abstract type Object end
abstract type Morphism end
```

## Categories

Categories without any additional structure do not need any 
fields or methods. We follow the example of the category of finite sets.

```@example Ex
struct FinSets <: Category end
```

## Objects

Objects need at least one field `parent` for the parent category or a method `parent`returning the respective category. Any other
information needed to work with the objects is arbitrary.

```@example Ex
struct FinSetObject <: Object
    parent::FinSets
    set::Set
end
```

Here we wrap a set to an object.

## Morphisms

Morphisms need to provide fields

- `domain`
- `codomain`

or methods

- `domain`
- `codomain`

For the category of sets we get

```@example Ex
struct FinSetMorphism <: Morphism
    domain::FinSetObject
    codomain::FinSetObject
    map
end
```

Where one now can design any coding for a morphism of sets
that fit the desired purpose.

## Required Methods

Necessary methods to implemented for morphisms, objects and categories are

- `compose(f::YourMorphism, g::YourMorphism)::YourMorphism` returning the composition ``g \circ f``.
- `id(X::YourObject)::YourMorphism` returning the identity morphism on ``X``.
- `Hom(X::YourObject, Y::YourObject)::AbstractHomSet` constructing an object `<:AbstractHomSet`.

Here anything extending `AbstractHomspace` needs to provide the fields 

- `domain`
- `codomain`

or methods

- `domain`
- `codomain`.

## Additional methods

Your category might have more structure. For categories which are 

- [abelian](AbelianCategories.md)
- [monoidal](MonoidalCategories.md)
- [tensor](TensorCategories.md)

visit the corresponding chapter. The interface supports additionally the following constructions and operations.

### (Co)Products

For a list of objects $X_1,...,X_n$ methods for the product shall
return an object representing the categorical product $\prod X_i$ together with the projection maps $p_i \colon \prod X_i \to X_i$. 

```julia
product(X::YourObject...)::Tuple{YourObject, Vector{YourMorphism}}
```

You may only implement a binary version `product(X,Y)` in which case TensorCategories extends it automatically to a list-version. Keep in mind that this might be devastating to the runtime, since iteratively applying a binary product involves composing morphisms which most likely is expensive.

TensorCategories will also generate an infix operator 

```julia 
×(X::YourObject, Y::YourObject)::YourObject
```
that only returns the object in question.

Dually the same applies for a coproduct $∐ Xᵢ$.

```julia
coproduct(X::YourObject...)::Tuple{YourObject, Vector{YourMorphism}}
```

### Initial and Terminal Objects

If a category has an initial and/or terminal object one can provide those.

```julia
initial_object(::Category)::Object
terminal_object(::Category)::Object
```


