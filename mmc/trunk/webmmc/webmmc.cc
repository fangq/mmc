/// Copyright (c) 2012 The Native Client Authors. All rights reserved.
/// Use of this source code is governed by a BSD-style license that can be
/// found in the LICENSE file.
///
/// @file webmmc.cc
/// This example demonstrates loading, running and scripting a very simple NaCl
/// module.  To load the NaCl module, the browser first looks for the
/// CreateModule() factory method (at the end of this file).  It calls
/// CreateModule() once to load the module code from your .nexe.  After the
/// .nexe code is loaded, CreateModule() is not called again.

#include <cstdio>
#include <string>
#include "ppapi/cpp/instance.h"
#include "ppapi/cpp/module.h"
#include "ppapi/cpp/var.h"
#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"
#include "tictoc.h"

class WebMMCObject : public pp::Instance {
 public:
  mcconfig cfg;
  tetmesh mesh;
  raytracer tracer;
  visitor master;
  double Eabsorb;
  RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  unsigned int i;
  float raytri;
  unsigned int threadid,ncomplete,t0,dt;

  explicit WebMMCObject(PP_Instance instance) : pp::Instance(instance){
      WebMMC_Reset();
  }
  virtual ~WebMMCObject() {
      mcx_clearcfg(&cfg);
  }
  int WebMMC_Reset(){
      //memset(master,0,sizeof(visitor));
      raytri=0.f;
      ncomplete=0;
      Eabsorb=0.0;
      mcx_clearcfg(&cfg);
      mcx_initcfg(&cfg);
      return 0;
  }

  virtual void HandleMessage(const pp::Var& var_message) {
    std::string msg = var_message.AsString();
    if(msg.find("MMC_MSG_RESET")==0){
         WebMMC_Reset();
    }else if(msg.find("{")==0){
         cJSON *jroot = cJSON_Parse(msg.c_str());
         if(jroot){
                mcx_loadjson(jroot,&cfg);
                cJSON_Delete(jroot);
         }
    }
  }
};

class WebMMCModule : public pp::Module {
 public:
  WebMMCModule() : pp::Module() {}
  virtual ~WebMMCModule() {}
  virtual pp::Instance* CreateInstance(PP_Instance instance) {
    return new WebMMCObject(instance);
  }
};

namespace pp {
Module* CreateModule() {
  return new WebMMCModule();
}
}  // namespace pp
