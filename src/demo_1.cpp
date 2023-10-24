#include <topology.h>

int main( int argc, char** argv )
{
    auto F = syc::topology::model::v_01::sketchable_forest<>(11);
    add_edge( 1, 6, F );
    add_edge( 6, 7, F );
    add_edge( 6, 8, F );
    add_edge( 3, 9, F );
    add_edge( 3, 10, F );
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

    return 0;
}
