#include "g4hrsRun.hh"
#include "g4hrsRunData.hh"

g4hrsRun *g4hrsRun::gSingleton = NULL;

g4hrsRun::g4hrsRun() {
    gSingleton = this;
    fRunData = new g4hrsRunData();
    fRunData->Init();
}

g4hrsRun::~g4hrsRun() {
}

g4hrsRun *g4hrsRun::GetRun() {
    if( gSingleton == NULL ) {
        gSingleton = new g4hrsRun();
    }
    return gSingleton;
}
