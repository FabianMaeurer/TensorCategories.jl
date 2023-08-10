# The Motivation 

This package began its journey asking the question "Can we play around 
with explicit categorical entities in the computer?".

By nature categorical operations and constructions are very generic and 
can be applied as long as the objects (or morphisms) are fitting. 
TensorCategories.jl provides an interface for categories with additional 
structure, precisely additive, linear, abelian, monoidal, tensor and 
fusion categories.

# Realizing Categories in The Computer

Due to the nature of category theory the realization of certain categories 
is very dependent on themselves. Thus the internal workings are generally 
up to the developer. As long as the interface for the desired additional
structures is implemented. 

Some kind of categories, i.e. fusion categories, are entirely described
(up to equivalence) by discrete data known as ``6j$-symbols``. Thus 
for such categories we can provide a datatype [`SixJCategory`](../SixJCategories/SixJCategories.md) 
to quickly work with categories given by such data.

# Mathematical Foundation

Throughout the package we will consider definitions and terminology as
provided in [EGNO](@cite).