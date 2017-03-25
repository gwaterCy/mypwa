#include "PWA_CTRL.H"

using namespace std;

void string_to_vector(string str, vector<string> &vec)
{
    while(str.find(',') != string::npos) {
        int i  = str.find(',');
        vec.push_back(str.substr(0, i));
        str = str.substr(i + 1);
    }
    vec.push_back(str);
}

static bool readConfigFile(const string & cfgfilepath, const string & key, string & value)
{
    fstream cfgFile;
    //cfgFile.open("dfa");//打开文件
    const char* fname = cfgfilepath.c_str();
    cfgFile.open(fname);//打开文件
    if( ! cfgFile.is_open())
    {
        cout<<"can not open cfg file!"<<endl;
        return false;
    }
    char tmp[1000];
    cout << "Begin readConfigFile" << endl;
    while(!cfgFile.eof())//循环读取每一行
    {
        cfgFile.getline(tmp,1000);//每行读取前1000个字符，前
        string line(tmp);
        //cout << "LINE:" << __LINE__ << tmp << endl;
        size_t pos = line.find('=');
        if(pos==string::npos) continue;
        string tmpKey = line.substr(0, pos);
        if(key==tmpKey)
        {
            value = line.substr(pos+1);//取=号之后
            cout << "value = " << value << endl;
            return true;
        }
    }
    return false;
}
void refresh_pwa_ctrl(const string config_file, PWA_CTRL & pwa_ctrl)
{
    string paraConstStr;
    readConfigFile(config_file, "const_paras", paraConstStr);
    string_to_vector(paraConstStr, pwa_ctrl.paraConstList);

    string actResStr;
    readConfigFile(config_file, "active_resonances", actResStr);
    string_to_vector(actResStr,pwa_ctrl.actResList);
    string actResStrKK;
    readConfigFile(config_file, "active_resonances_kk", actResStrKK);
    if (actResStrKK != "NONE") {
        string_to_vector(actResStrKK,pwa_ctrl.actResListKK);
    }
    string actResStrPP;
    readConfigFile(config_file, "active_resonances_pp", actResStrPP);
    if (actResStrPP != "NONE") {
        string_to_vector(actResStrPP,pwa_ctrl.actResListPP);
    }

    readConfigFile(config_file, "phipp_phsp_file", pwa_ctrl.phsp_phipp);
    readConfigFile(config_file, "phikk_phsp_file", pwa_ctrl.phsp_phikk);
    readConfigFile(config_file, "phipp_data_file", pwa_ctrl.data_phipp);
    readConfigFile(config_file, "phikk_data_file", pwa_ctrl.data_phikk);

    readConfigFile(config_file, "reload_paras_file", pwa_ctrl.reloadFitParaFile);

    readConfigFile(config_file, "work_path", pwa_ctrl.workPath);
    readConfigFile(config_file, "store_paras_file", pwa_ctrl.storeFitParaFile);
    pwa_ctrl.storeFitParaFile = pwa_ctrl.workPath + pwa_ctrl.storeFitParaFile;
    readConfigFile(config_file, "outf_phipp_file", pwa_ctrl.outf_phipp);
    pwa_ctrl.outf_phipp = pwa_ctrl.workPath + pwa_ctrl.outf_phipp;
    readConfigFile(config_file, "outf_phikk_file", pwa_ctrl.outf_phikk);
    pwa_ctrl.outf_phikk = pwa_ctrl.workPath + pwa_ctrl.outf_phikk;
    readConfigFile(config_file, "proj_phipp_file", pwa_ctrl.proj_phipp);
    pwa_ctrl.proj_phipp = pwa_ctrl.workPath + pwa_ctrl.proj_phipp;
    readConfigFile(config_file, "proj_phikk_file", pwa_ctrl.proj_phikk);
    pwa_ctrl.proj_phikk = pwa_ctrl.workPath + pwa_ctrl.proj_phikk;

    readConfigFile(config_file, "phipp_index_weight_file", pwa_ctrl.idx_pp);
    readConfigFile(config_file, "phikk_index_weight_file", pwa_ctrl.idx_kk);

    cout << "phsp_phipp = " << pwa_ctrl.phsp_phipp << endl;
    cout << "data_phipp = " << pwa_ctrl.data_phipp << endl;
    cout << "idx_pp = " << pwa_ctrl.idx_pp << endl;
}

