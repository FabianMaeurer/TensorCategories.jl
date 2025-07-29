

@testset "AnyonWiki" begin 

    # Test the construction 
    @testset "Construction Categories" begin

        # Test loading of simple categories
        @test try 
            C = anyonwiki(3,1,2,1,3,0,1)
            true
        catch
            false
        end
        @test try 
            C = anyonwiki(4,1,2,2,4,0,2)
            true
        catch
            false
        end
    end

    # Test center loading 
    @testset "Centers of anyonwiki" begin

        # Test loading of simple centers
        @test try 
            C = anyonwiki_center(3,1,2,1,3,0,1)
            true
        catch
            false
        end
        @test try 
            C = anyonwiki_center(4,1,2,2,4,0,2)
            true
        catch
            false
        end
    end

end
