

@testset "AnyonWiki" begin 

    # Test the construction 
    @testset "Construction Categories" begin

        # Test loading of simple categories
        @test try 
            C = anyonwiki(2,1,2,1,3,0,1)
            true
        catch
            false
        end
    end

    # Test center loading 
    @testset "Centers of anyonwiki" begin

        # Test loading of simple centers
        @test try 
            C = anyonwiki_center(2,1,2,1,3,0,1)
            true
        catch
            false
        end
    end

end
