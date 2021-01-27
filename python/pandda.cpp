// Copyright 2018 Global Phasing Ltd.

#include <iostream>

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped
#include "gemmi/neighbor.hpp"
#include "gemmi/tostr.hpp"
#include "gemmi/fourier.hpp"  // for get_f_phi_on_grid

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "common.h"  // for normalize_index

#include "ThreadPool.h"


namespace py = pybind11;
using namespace gemmi;


std::vector<double> interpolate_to_list(
    Grid<float> moving_map,
    Grid<float> interpolated_map, 
    std::vector<std::vector<int>> point_vec,
    std::vector<std::vector<double>> pos_vec,
    std::vector<Transform> transform_vec,
    std::vector<std::vector<double>> com_moving_vec,
    std::vector<std::vector<double>> com_reference_vec
)
{

    std::vector<double> vals(point_vec.size());

    for (std::size_t i=0; i < point_vec.size(); i++)
    {
        // Position
        std::vector<int> point = point_vec[i];

        Fractional fractional = Fractional(
            point[0] * (1.0 / interpolated_map.nu), 
            point[1] * (1.0 / interpolated_map.nv), 
            point[2] * (1.0 / interpolated_map.nw)
            );
        Position pos = interpolated_map.unit_cell.orthogonalize(fractional);
        Transform transform = transform_vec[i];
        std::vector<double> com_moving = com_moving_vec[i];
        std::vector<double> com_reference = com_reference_vec[i];

        //Subtract reference com
        pos.x -= com_reference[0];
        pos.y -= com_reference[1];
        pos.z -= com_reference[2];

        //transform
        Position pos_moving = Position(transform.apply(pos));

        // add moving com
        pos_moving.x += com_moving[0];
        pos_moving.y += com_moving[1];
        pos_moving.z += com_moving[2];

        // fractionalise
        Fractional pos_moving_fractional = moving_map.unit_cell.fractionalize(pos_moving);

        // interpolate
        float interpolated_value = moving_map.interpolate_value(pos_moving_fractional);

        // assign
        vals[i] = interpolated_value;


    };

    return vals;


}

Grid<float> interpolate_points(
    Grid<float> moving_map,
    Grid<float> interpolated_map, 
    std::vector<std::vector<int>> point_vec,
    std::vector<std::vector<double>> pos_vec,
    std::vector<Transform> transform_vec,
    std::vector<std::vector<double>> com_moving_vec,
    std::vector<std::vector<double>> com_reference_vec
    )
{
    for (std::size_t i=0; i < point_vec.size(); i++)
    {

        // Release gil for threading support
        py::gil_scoped_release release;

        // Position
        std::vector<int> point = point_vec[i];
        std::cout << "Point: " << point[0] << " " << point[1] << " " << point[2] << "\n";

        Fractional fractional = Fractional(
            point[0] * (1.0 / interpolated_map.nu), 
            point[1] * (1.0 / interpolated_map.nv), 
            point[2] * (1.0 / interpolated_map.nw)
            );
        Position pos = interpolated_map.unit_cell.orthogonalize(fractional);
        Transform transform = transform_vec[i];
        std::cout << "Trasform: " << transform.vec.x << " " << transform.vec.y << " " << transform.vec.z << "\n";
        std::cout << "Trasform: " << transform.mat[0][0] << " " << transform.mat[0][1] << " " << transform.mat[0][2] << "\n";
        std::cout << "Trasform: " << transform.mat[1][0] << " " << transform.mat[1][1] << " " << transform.mat[1][2] << "\n";
        std::cout << "Trasform: " << transform.mat[2][0] << " " << transform.mat[2][1] << " " << transform.mat[2][2] << "\n";
        

        std::vector<double> com_moving = com_moving_vec[i];
        std::cout << "com_moving: " << com_moving[0] << " " << com_moving[1] << " " << com_moving[2] << "\n";
        std::vector<double> com_reference = com_reference_vec[i];
        std::cout << "com_reference: " << com_reference[0] << " " << com_reference[1] << " " << com_reference[2] << "\n";

        //Subtract reference com
        pos.x -= com_reference[0];
        pos.y -= com_reference[1];
        pos.z -= com_reference[2];

        //transform
        Position pos_moving = Position(transform.apply(pos));

        // add moving com
        pos_moving.x += com_moving[0];
        pos_moving.y += com_moving[1];
        pos_moving.z += com_moving[2];

        // fractionalise
        Fractional pos_moving_fractional = moving_map.unit_cell.fractionalize(pos_moving);

        // interpolate
        float interpolated_value = moving_map.interpolate_value(pos_moving_fractional);

        // assign
        interpolated_map.set_value(
            point[0],
            point[1],
            point[2],
            interpolated_value
            );


    };

    return interpolated_map;


}


// std::vector<Fractional> points_to_positions(
//     Grid<float> interpolated_map, 
//     std::vector<std::vector<int>> point_vec
//     )
// {

//     std::vector<Fractional> position_vec(position_vec.size());

//     for (std::size_t i=0; i < point_vec.size(); i++)
//     {
//         // Position
//         std::vector<int> point = point_vec[i];

//         Fractional fractional = Fractional(
//             point[0] * (1.0 / interpolated_map.nu), 
//             point[1] * (1.0 / interpolated_map.nv), 
//             point[2] * (1.0 / interpolated_map.nw)
//             );
//         Position pos = interpolated_map.unit_cell.orthogonalize(fractional);
        
//         position_vec[i] = pos;

//     }

//     return position_vec;

// }

// std::vector<Fractional> transform_positions(
//     position_vec,
//     transform_vec,
//     com_moving_vec,
//     com_reference_vec
// )
// {
//     std::vector<Fractional> fractional_vec(position_vec.size());
    
//     for (std::size_t i=0; i < position_vec.size(); i++)
//     {
//         // Position
//         Position pos = position_vec[i];
//         Transform transform = transform_vec[i];
//         std::vector<double> com_moving = com_moving_vec[i];
//         std::vector<double> com_reference = com_reference_vec[i];

//         //Subtract reference com
//         pos.x -= com_reference[0];
//         pos.y -= com_reference[1];
//         pos.z -= com_reference[2];

//         //transform
//         Position pos_moving = Position(transform.apply(pos));

//         // add moving com
//         pos_moving.x += com_moving[0];
//         pos_moving.y += com_moving[1];
//         pos_moving.z += com_moving[2];

//         // fractionalise
//         Fractional pos_moving_fractional = moving_map.unit_cell.fractionalize(pos_moving);

//         fractional_vec[i] = pos_moving_fractional;

//     }

//     return fractional_vec;
// }

// Grid<float> interpolate_points_from_positions(
//     Grid<float> moving_map,
//     Grid<float> interpolated_map, 
//     std::vector<std::vector<int>> point_vec,
//     std::vector<Fractional> position_vec
//     )
// {
//     // Release gil for threading support
//     py::gil_scoped_release release;

//     // Iterate over the corresponding points and positions to interpolate at 
//     for (std::size_t i=0; i < position_vec.size(); i++)
//     {
//         std::vector<int> point = point_vec[i];

//         Fractional position = position_vec[i];

//         // interpolate
//         float interpolated_value = moving_map.interpolate_value(position);

//         // assign
//         interpolated_map.set_value(
//             point[0],
//             point[1],
//             point[2],
//             interpolated_value
//             );


//     };

//     return interpolated_map;

// }



// Grid<float> interpolate_points(
//     Grid<float> moving_map,
//     Grid<float> interpolated_map, 
//     std::vector<std::vector<int>> point_vec,
//     std::vector<std::vector<double>> pos_vec,
//     std::vector<Transform> transform_vec,
//     std::vector<std::vector<double>> com_moving_vec,
//     std::vector<std::vector<double>> com_reference_vec
//     )
// {



//     Grid<float> interpolated_map = interpolate_points_from_positions(
//             Grid<float> moving_map,
//             Grid<float> interpolated_map, 
//             std::vector<std::vector<int>> point_vec,
//             std::vector<Fractional> position_vec
//             );

//     return interpolated_map;
// }

// float get_rmsd(
//     std::vector<float>,
//     std::vector<float>,
// )
// {


// }

// float optimise_stochastic_downhill(
//     std::function<float(float)> func,
//     float step,
//     float start
// )
// {

// }

// std::vector<float> apply_exponential_scale(
//     std::vector<float> array,
// )
// {


// }

// struct Solver
// {


//     std::vector<float> reference_array;
//     std::vector<float> moving_array;
//     std::vector<float> working_array;
//     std::vector<float> resolution_array;


//     void operator()(float scale)
//         {
            
//             reset_working_array(working_array, moving_array);

//             apply_exponential_scale(scale, working_array);

//             float rmsd = get_rmsd(reference_array);

//             return rmsd; 

//         }
// };

// float calculate_isotropic_b_factor_scale(
//     std::vector<float> reference_array,
//     std::vector<float> moving_array,
//     std::vector<float> resolution_array
// )
// {

//     Solver solver;
//     solver.reference_array ;
//     solver.moving_array ;
//     solver.working_array ;
//     solver.resolution_array ;
    

//     scaling = optimise_stochastic_downhill(solver, 0.05, 0.0);

//     scaled_array = apply_exponential_scale(scaling, moving_array);

//     return scaling; 

// }


void add_pandda(py::module& m) {
    m.def(
        "interpolate_points", 
        &interpolate_points, 
        "Interpolates a list of points."
    );

    m.def(
        "interpolate_to_list", 
        &interpolate_to_list, 
        "Interpolates a list of points."
    );

}

