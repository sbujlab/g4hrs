#ifndef __REMOLLRUN_HH
#define __REMOLLRUN_HH

/*!
 * All the information on the run
 * The data object will get put into the output
 * stream

   This is implemented in the soliton model
 */

#include "g4hrsRunData.hh"

class g4hrsRun {

private:
    static g4hrsRun *gSingleton;
    g4hrsRun();

    g4hrsRunData *fRunData;

public:
    static g4hrsRun *GetRun();
    ~g4hrsRun();

    g4hrsRunData *GetData() {
        return fRunData;
    }
};

#endif//__REMOLLRUN_HH
