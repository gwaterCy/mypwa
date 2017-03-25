#include "../codes_frac/fitproxy.h"

int main()
{
    fitproxy *A = new fitproxy();
    A->init_all_files_name();
    A->init_input_argset();
    A->setup_resonances();
    // A->add_res0_list("f01750", 1.7, 2.2 0.15, 0.7) here!!!
    A->createlist_allparas();
    A->gp("f01000a_mss_")->setConstant();
    A->gp("f01000a_wdt_")->setConstant();
    A->gp("f21000a_mss_")->setConstant();
    A->gp("f21000a_wdt_")->setConstant();
    A->print_all_paras();

    A->read_data();
    A->init_pdf();
    A->act_res980("f00980");
    A->act_res0("f01000");
    A->act_res0("f01500");
//  A->act_res2("f21270");
    A->act_res2("f21525");
    //A->act_res1m("1m1800");
    //A->act_res1p("1p1800");

    //A->act_res0("f01500");
    //    A->act_res2("f21270");
    //    A->act_res2("f21525");
    // A->act_res0("f01750");  here!!!
    A->print_fit_paras();

    A->FIT();
    cout << "+++++print_fit_paras" << endl;
    A->print_fit_paras();
    cout << "++++++store_fit_paras" << endl;
    A->store_fit_paras("fit_base.root");
    cout << "Fit is done!" << endl;
    delete A;
}


