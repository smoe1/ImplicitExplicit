#ifndef _VPDATACART2_H_
#define _VPDATACART2_H_

class dTensor2;
class dTensorBC3;
class dTensorBC4;

class VPDataCart2
{
    public:
        dTensor2*   node1d;
        dTensorBC3* Efield;
        dTensorBC3* Efield_old;
        dTensorBC3* v1d;
        dTensorBC3* KE;
        dTensorBC3* KE_old;
        dTensorBC4* qold;

        // Default constructor
        VPDataCart2():is_initialized(false){}
        void init();

        // Destructor
        ~VPDataCart2();

    private:
        
        // Flag to determine whether or not class has been initialized
        bool is_initialized;

    // TODO - where is the destructor?  This "class" will have memory leaks if
    // instantiated multiple times.

};
extern VPDataCart2 vpDataCart2;

#endif

