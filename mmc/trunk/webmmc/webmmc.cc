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

class WebMMCInstance : public pp::Instance {
 public:
  explicit WebMMCInstance(PP_Instance instance) : pp::Instance(instance)
  {}
  virtual ~WebMMCInstance() {}

  virtual void HandleMessage(const pp::Var& var_message) {
    //std::string msg = var_message.AsString();
    //char *jsonmsg=msg.c_str();
  }
};

class WebMMCModule : public pp::Module {
 public:
  WebMMCModule() : pp::Module() {}
  virtual ~WebMMCModule() {}
  virtual pp::Instance* CreateInstance(PP_Instance instance) {
    return new WebMMCInstance(instance);
  }
};

namespace pp {
Module* CreateModule() {
  return new WebMMCModule();
}
}  // namespace pp
