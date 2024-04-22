# Generic Methods 

Many constructions for additive, abelian, linear or monoidal categories are entangled. Thus we provide a 
vast list of methods to compute objects or morphisms using other methods. 

Keep in mind that the performance will usually be much better if the following methods are 
overwritten form custom types. 

```@autodocs
Modules = [TensorCategories]
Pages = ["AbstractMethods.jl", "Fallbacks.jl"]
Private = false
```