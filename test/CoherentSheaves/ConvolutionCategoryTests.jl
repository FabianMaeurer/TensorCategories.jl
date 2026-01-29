
@testset "ConvolutionCategory Tests" begin
#     # Test ConvolutionCategory constructor
#     @testset "Constructor Tests" begin
#         G = symmetric_group(3)
#         K,_ = finite_field(2)
#         X = gset(G, [1, 2, 3])
#         conv_cat = convolution_category(K, X)
        
#         @test conv_cat.group == G
#         @test conv_cat.base_ring == K
#     end

#     # Test TensorCategories.ConvolutionObject constructor
#     @testset "Object Tests" begin
#         G = symmetric_group(3)
#         K,_ = finite_field(2)
#         X = gset(G, [1, 2, 3])
#         conv_cat = convolution_category(K, X)
#         sheaf = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [zero(representation_category(K, G))])
#         conv_obj = TensorCategories.ConvolutionObject(sheaf, conv_cat)
        
#         @test conv_obj.sheaf == sheaf
#         @test conv_obj.parent == conv_cat
#     end

#     # Test TensorCategories.ConvolutionMorphism constructor
#     @testset "Morphism Tests" begin
#         G = symmetric_group(3)
#         K,_ = finite_field(2)
#         X = gset(G, [1, 2, 3])
#         conv_cat = convolution_category(K, X)
#         sheaf = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [zero(representation_category(K, G))])
#         conv_obj = TensorCategories.ConvolutionObject(sheaf, conv_cat)
#         morphism = TensorCategories.CohSheafMorphism(sheaf, sheaf, zero_matrix(K, 1, 1))
#         conv_morph = TensorCategories.ConvolutionMorphism(conv_obj, conv_obj, morphism)
        
#         @test conv_morph.domain == conv_obj
#         @test conv_morph.codomain == conv_obj
#         @test conv_morph.m == morphism
#     end

    # Test is_multitensor and is_fusion
    @testset "Category Properties Tests" begin
        G = symmetric_group(3)
        K,_ = finite_field(2)
        X = gset(G, [1, 2, 3])
        conv_cat = convolution_category(K, X)
        
        @test is_multitensor(conv_cat)
        @test is_fusion(conv_cat) == (mod(order(G), characteristic(K)) != 0)
    end

    # # Test tensor product of ConvolutionObjects
    # @testset "Tensor Product Tests" begin
    #     G = symmetric_group(3)
    #     K,_ = finite_field(2)
    #     X = gset(G, [1, 2, 3])
    #     conv_cat = convolution_category(K, X)
    #     sheaf1 = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [zero(representation_category(K, G))])
    #     sheaf2 = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [one(representation_category(K, G))])
    #     conv_obj1 = TensorCategories.ConvolutionObject(sheaf1, conv_cat)
    #     conv_obj2 = TensorCategories.ConvolutionObject(sheaf2, conv_cat)
        
    #     tensor_prod = tensor_product(conv_obj1, conv_obj2)
        
    #     @test tensor_prod.parent == conv_cat
    #     @test tensor_prod.sheaf == conv_cat.projectors[2](conv_cat.projectors[1](sheaf1) ⊗ conv_cat.projectors[3](sheaf2))
    # end

    # # Test tensor product of ConvolutionMorphisms
    # @testset "Tensor Product Morphism Tests" begin
    #     G = symmetric_group(3)
    #     K,_ = finite_field(2)
    #     X = gset(G, [1, 2, 3])
    #     conv_cat = convolution_category(K, X)
    #     sheaf1 = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [zero(representation_category(K, G))])
    #     sheaf2 = TensorCategories.CohSheafObject(conv_cat.squaredCoh, [one(representation_category(K, G))])
    #     conv_obj1 = TensorCategories.ConvolutionObject(sheaf1, conv_cat)
    #     conv_obj2 = TensorCategories.ConvolutionObject(sheaf2, conv_cat)
    #     morphism1 = TensorCategories.CohSheafMorphism(sheaf1, sheaf1, zero_matrix(K, 1, 1))
    #     morphism2 = TensorCategories.CohSheafMorphism(sheaf2, sheaf2, one_matrix(K, 1, 1))
    #     conv_morph1 = TensorCategories.ConvolutionMorphism(conv_obj1, conv_obj1, morphism1)
    #     conv_morph2 = TensorCategories.ConvolutionMorphism(conv_obj2, conv_obj2, morphism2)
        
    #     tensor_prod_morph = tensor_product(conv_morph1, conv_morph2)
        
    #     @test tensor_prod_morph.domain == tensor_product(conv_obj1, conv_obj2)
    #     @test tensor_prod_morph.codomain == tensor_product(conv_obj1, conv_obj2)
    #     @test tensor_prod_morph.m == conv_cat.projectors[2](conv_cat.projectors[1](morphism1) ⊗ conv_cat.projectors[3](morphism2))
    # end

    # Test simple objects
    @testset "Simple Objects Tests" begin
        G = symmetric_group(3)
        K,_ = finite_field(2)
        X = gset(G, [1, 2, 3])
        conv_cat = convolution_category(K, X)
        
        simpls = simples(conv_cat)
        @test length(simpls) > 0
        @test all(is_simple, simpls)
        
        for simple_obj in simpls
            @test is_simple(simple_obj)
        end
    end

end
