#ifndef _POISSONPARAMS_H_
#define _POISSONPARAMS_H_

class IniDocument;
class PoissonParams
{  
    public:
        int lft_bc_type;
        int rgt_bc_type;
        double lft_val;
        double rgt_val;
        void init(IniDocument& ini_doc);
};
extern PoissonParams poissonParams;

#endif
