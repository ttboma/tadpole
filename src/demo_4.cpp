#include <fstream>
#include <topology.h>

int main( int argc, char** argv )
{
    auto F = syc::topology::model::v_01::sketchable_forest<>(9);
    add_edge( 3, 7, F );
    add_edge( 7, 8, F );
    add_edge( 8, 9, F );

    F.make_sketch();
    std::cout << "current heand: none\n";
    F.txt_vitualizer();

    while (true)
    {
        std::cout << "please enter two **valid** index to make slice (enter q to exit):\n";
        if ( int u, v, h; std::cin >> u >> v )
        {
            h = F.make_slice( u, v );
            std::cout << "current heand: " << h << "\n";
            F.txt_vitualizer();
        }
        else 
        {
            break;
        }
    }
    
    auto out = std::ofstream("demo3.txt");
    F.txt_vitualizer2( out );

    return 0;
}
