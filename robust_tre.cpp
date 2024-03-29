#include <ppl.hh>
#include <gmpxx.h>
#include <pybind11/pybind11.h>

#include "bound.hpp"
#include "zone.hpp"
#include "zone_set.hpp"
#include "utils.hpp"

using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;

int add(int i, int j) {
    return i + j;
}

/* Fully accurate when zones don't intersect */
/* Gives a conservative estimate otherwise */
template <typename T>
timedrel::zone_set<T> time_robust_match_translation(timedrel::zone_set<T> &zs_in, T r_lbound){
    timedrel::zone_set<T> zs_res;

    Variable x(0),y(1),delta(2);

    /* Convert zones to robustness polyhedra */
    for(auto z : zs_in){
        Constraint_System cs;
        mpq_class x_min(z.get_bmin().value);
        mpq_class x_max(z.get_bmax().value);
        mpq_class y_min(z.get_emin().value);
        mpq_class y_max(z.get_emax().value);
        mpq_class delta_min(z.get_dmin().value);
        mpq_class delta_max(z.get_dmax().value);
        mpq_class delta_lbound(r_lbound);

        /* Add zone constraints */
        cs.insert(x_min.get_num() <= x*x_min.get_den());
        cs.insert(y_min.get_num() <= y*y_min.get_den());
        cs.insert(x_max.get_den()*x <= x_max.get_num());
        cs.insert(y_max.get_den()*y <= y_max.get_num());
        cs.insert(delta_min.get_num() <= y*delta_min.get_den()-x*delta_min.get_den());
        cs.insert(y*delta_max.get_den()-x*delta_max.get_den() <= delta_max.get_num());

        /* Add robustness constraints */
        cs.insert(delta >= 0);
        cs.insert(delta*x_min.get_den() <= x*x_min.get_den()-x_min.get_num());
        cs.insert(delta*y_min.get_den() <= y*y_min.get_den()-y_min.get_num());
        cs.insert(delta*x_max.get_den() <= x_max.get_num()-x*x_max.get_den());
        cs.insert(delta*y_max.get_den() <= y_max.get_num()-y*y_max.get_den());
        cs.insert(delta*delta_lbound.get_den() >= delta_lbound.get_num());

        C_Polyhedron phedra(cs);
        phedra.unconstrain(delta);

        Coefficient rx_min_n, rx_min_d;
        Coefficient ry_min_n, ry_min_d;
        Coefficient rd_min_n, rd_min_d;

        Coefficient rx_max_n, rx_max_d;
        Coefficient ry_max_n, ry_max_d;
        Coefficient rd_max_n, rd_max_d;

        bool minim, maxim;
        phedra.minimize(x, rx_min_n, rx_min_d, minim);
        phedra.minimize(y, ry_min_n, ry_min_d, minim);
        phedra.minimize(y-x, rd_min_n, rd_min_d, minim);

        phedra.maximize(x, rx_max_n, rx_max_d, maxim);
        phedra.maximize(y, ry_max_n, ry_max_d, maxim);
        phedra.maximize(y-x, rd_max_n, rd_max_d, maxim);

        mpq_class rx_min_q(rx_min_n, rx_min_d);
        mpq_class ry_min_q(ry_min_n, ry_min_d);
        mpq_class rd_min_q(rd_min_n, rd_min_d);
        mpq_class rx_max_q(rx_max_n, rx_max_d);
        mpq_class ry_max_q(ry_max_n, ry_max_d);
        mpq_class rd_max_q(rd_max_n, rd_max_d);

        T rx_min = rx_min_q.get_d();
        T ry_min = ry_min_q.get_d();
        T rd_min = rd_min_q.get_d();
        T rx_max = rx_max_q.get_d();
        T ry_max = ry_max_q.get_d();
        T rd_max = rd_max_q.get_d();

        zs_res.add({rx_min, rx_max, ry_min, ry_max, rd_min, rd_max}, {1,1,1,1,1,1});
    }

    return zs_res;
}

/* Fully accurate when zones don't intersect */
/* Gives a conservative estimate otherwise */
template <typename T>
T get_time_robustness_translation(timedrel::zone_set<T> &zs_in, T l, T u){
    Variable x(0),y(1),delta(2);

    T rob_value = 0;

    /* Convert zones to robustness polyhedra */
    for(auto z : zs_in){
        Constraint_System cs;
        mpq_class x_min(z.get_bmin().value);
        mpq_class x_max(z.get_bmax().value);
        mpq_class y_min(z.get_emin().value);
        mpq_class y_max(z.get_emax().value);
        mpq_class delta_min(z.get_dmin().value);
        mpq_class delta_max(z.get_dmax().value);
        mpq_class px(l);
        mpq_class py(u);

        /* Add zone constraints */
        cs.insert(x_min.get_num() <= x*x_min.get_den());
        cs.insert(y_min.get_num() <= y*y_min.get_den());
        cs.insert(x_max.get_den()*x <= x_max.get_num());
        cs.insert(y_max.get_den()*y <= y_max.get_num());
        cs.insert(delta_min.get_num() <= y*delta_min.get_den()-x*delta_min.get_den());
        cs.insert(y*delta_max.get_den()-x*delta_max.get_den() <= delta_max.get_num());

        /* Add robustness constraints */
        cs.insert(delta >= 0);
        cs.insert(delta*x_min.get_den() <= x*x_min.get_den()-x_min.get_num());
        cs.insert(delta*y_min.get_den() <= y*y_min.get_den()-y_min.get_num());
        cs.insert(delta*x_max.get_den() <= x_max.get_num()-x*x_max.get_den());
        cs.insert(delta*y_max.get_den() <= y_max.get_num()-y*y_max.get_den());

        /* Add interval time point constraints */
        cs.insert(x*px.get_den() == px.get_num());
        cs.insert(y*py.get_den() == py.get_num());

        C_Polyhedron phedra(cs);

        if(!phedra.is_empty()){
            Coefficient r_max_n, r_max_d;
            bool maxim;
            phedra.maximize(delta, r_max_n, r_max_d, maxim);
            mpq_class rob_max(r_max_n, r_max_d);
            if(rob_max.get_d() > rob_value){
                rob_value = rob_max.get_d();
            }
        }
    }

    return rob_value;
}

/* Fully accurate (hopefully) */
template <typename T>
T get_time_robustness_translation_optimal(timedrel::zone_set<T> &zs_in, T l, T u, T scope_start, T scope_end){
    T rob_value_right = 0, rob_value_left = 0, rob_value = 0;
    auto zs_line = timedrel::zone_set<T>();
    zs_line.add({scope_start,scope_end,scope_start,scope_end,u-l,u-l},{1,1,1,1,1,1});

    auto zs_inter = timedrel::zone_set<T>::intersection(zs_in, zs_line);

    std::vector<T> border_points_right, border_points_left;
    std::vector<T> eborder_points_right, eborder_points_left;
    /* Get points to the right and points to the left */
    for(auto z : zs_inter){
        T sp = z.get_bmin().value;
        T ep = z.get_bmax().value;
        T esp = z.get_emin().value;
        T eep = z.get_emax().value;
        if(sp >= l){
            border_points_right.push_back(sp);
            eborder_points_right.push_back(esp);
        }
        if(ep >= l){
            border_points_right.push_back(ep);
            eborder_points_right.push_back(eep);
        }
        if(sp <= l){
            border_points_left.push_back(sp);
            eborder_points_left.push_back(esp);
        }
        if(ep <= l){
            border_points_left.push_back(ep);
            eborder_points_left.push_back(eep);
        }
    }
    sort(border_points_right.begin(), border_points_right.end());
    sort(border_points_left.begin(), border_points_left.end());

    T old_point = l;
    T new_point = l;
    T eold_point = u;
    T enew_point = u;
    /* Compute robustness to the right */
    for(int i=0; i < border_points_right.size(); i++){
        new_point = border_points_right[i];
        enew_point = eborder_points_right[i];
        auto zs_segment = timedrel::zone_set<T>();
        zs_segment.add({old_point, new_point,
            eold_point, enew_point,
            u-l,u-l},{1,1,1,1,1,1});

        if(timedrel::zone_set<T>::includes(zs_inter, zs_segment)){
            old_point = new_point;
            eold_point = enew_point;
        }else{
            break;
        }
    }
    /* Assign robustness value to the right */
    rob_value_right = old_point - l;

    old_point = l;
    new_point = l;
    eold_point = u;
    enew_point = u;
    /* Compute robustness to the left */
    int i = border_points_left.size() - 1;
    while(i >= 0){
        new_point = border_points_left[i];
        enew_point = eborder_points_left[i];
        auto zs_segment = timedrel::zone_set<T>();
        zs_segment.add({new_point, old_point,
            enew_point, eold_point,
            u-l,u-l},{1,1,1,1,1,1});

        if(timedrel::zone_set<T>::includes(zs_inter, zs_segment)){
            old_point = new_point;
            eold_point = enew_point;
        }else{
            break;
        }
        i--;
    }
    /* Assign robustness value to the left */
    rob_value_left = l - old_point;

    /* Assign robustness value */
    rob_value = (rob_value_right > rob_value_left)?rob_value_left:rob_value_right;

    return rob_value;
}

template <typename T>
void print_zone_set(timedrel::zone_set<T> &zs_in){
    std::cout<<"------"<<std::endl;
    for(auto z : zs_in){
        std::cout<<z<<std::endl;
    }
    std::cout<<"------"<<std::endl;
}

namespace py = pybind11;

using T = double;

PYBIND11_MODULE(robust_tre, m) {
    m.doc() = "timedrel robust plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");

    using namespace timedrel;

    typedef lower_bound<T> lower_bound_type;
    typedef upper_bound<T> upper_bound_type;

    m.def("trmtrans", &time_robust_match_translation<T>);
    m.def("zsetprint", &print_zone_set<T>);
    m.def("trobustness", &get_time_robustness_translation<T>);
    m.def("trobustness_opt", &get_time_robustness_translation_optimal<T>);

    py::class_<lower_bound_type>(m, "lower_bound")
        .def(py::init<T, bool>())
        .def_readonly("value", &lower_bound_type::value)
        .def_readonly("sign", &lower_bound_type::sign)
    ;

    m.def("lt", &upper_bound_type::strict);
    m.def("leq", &upper_bound_type::nonstrict);

    py::class_<upper_bound_type>(m, "upper_bound")
        .def(py::init<T, bool>())
        .def_readonly("value", &upper_bound_type::value)
        .def_readonly("sign", &upper_bound_type::sign)
    ;

    m.def("gt", &lower_bound_type::strict);
    m.def("geq", &lower_bound_type::nonstrict);

    typedef zone<T> zone_type;
    typedef zone_set<T> zone_set_type;

    py::class_<zone_type>(m, "zone")
        .def("bmin", &zone_type::get_bmin)
        .def("bmax", &zone_type::get_bmax)
        .def("emin", &zone_type::get_emin)
        .def("emax", &zone_type::get_emax)
        .def("dmin", &zone_type::get_dmin)
        .def("dmax", &zone_type::get_dmax)

        .def<zone_type (*)(
            const lower_bound_type&, const upper_bound_type&, 
            const lower_bound_type&, const upper_bound_type&, 
            const lower_bound_type&, const upper_bound_type&)>
        ("make", &zone_type::make)
    ;

    py::class_<zone_set_type>(m, "zone_set")
        .def(py::init<>())
        .def<void (zone_set_type::*)(const zone_type&)>("add", &zone_set_type::add)
        .def<void (zone_set_type::*)(const std::array<T, 6>&)>("add", &zone_set_type::add)
        .def<void (zone_set_type::*)(const std::array<T, 6>&, const std::array<bool, 6>&)>("add", &zone_set_type::add)
        .def("add_from_period", &zone_set_type::add_from_period)
        .def("add_from_period_rise_anchor", &zone_set_type::add_from_period_rise_anchor)
        .def("add_from_period_fall_anchor", &zone_set_type::add_from_period_fall_anchor)
        .def("add_from_period_both_anchor", &zone_set_type::add_from_period_both_anchor)
        .def("empty", &zone_set_type::empty)
        .def("__iter__", [](const zone_set_type &s) { return py::make_iterator(s.cbegin(), s.cend()); },
                         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    ;

    m.def("filter", &zone_set_type::filter);
    m.def("includes", &zone_set_type::includes);

    // Set operations
    m.def<zone_set_type (*)(const zone_set_type&)>("complementation", &zone_set_type::complementation);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("duration_restriction", &zone_set_type::duration_restriction);
    m.def<zone_set_type (*)(const zone_set_type&, const zone_set_type&)>("union", &zone_set_type::set_union);
    m.def<zone_set_type (*)(const zone_set_type&, const zone_set_type&)>("intersection", &zone_set_type::intersection);
    m.def<zone_set_type (*)(const zone_set_type&, const zone_set_type&)>("difference", &zone_set_type::set_difference);

    // Sequential operations
    m.def<zone_set_type (*)(const zone_set_type&, const zone_set_type&)>("concatenation", &zone_set_type::concatenation);
    m.def<zone_set_type (*)(const zone_set_type&)>("transitive_closure", &zone_set_type::transitive_closure);

    // Modal operations of the logic of time periods
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_starts", &zone_set_type::diamond_starts);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_started_by", &zone_set_type::diamond_started_by);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_finishes", &zone_set_type::diamond_finishes);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_finished_by", &zone_set_type::diamond_finished_by);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_meets", &zone_set_type::diamond_meets);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("diamond_met_by", &zone_set_type::diamond_met_by);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_starts", &zone_set_type::box_starts);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_started_by", &zone_set_type::box_started_by);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_finishes", &zone_set_type::box_finishes);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_finished_by", &zone_set_type::box_finished_by);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_meets", &zone_set_type::box_meets);
    m.def<zone_set_type (*)(const zone_set_type&, T, T)>("box_met_by", &zone_set_type::box_met_by);
}
