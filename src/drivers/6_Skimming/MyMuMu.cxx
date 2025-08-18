// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MyMuMu
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "MyAnalysis.h"
#include "CERNfunctions.h"
#include "vector"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace {
  void TriggerDictionaryInitialization_MyMuMu_Impl() {
    static const char* headers[] = {
"MyAnalysis.h",
"CERNfunctions.h",
0
    };
    static const char* includePaths[] = {
"/disk1/root/include/",
"/home/aliko/HmmAnalysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MyMuMu dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MyMuMu dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MyAnalysis.h"
#include "CERNfunctions.h"
#include "vector"
#ifdef __MAKECINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MyMuMu",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MyMuMu_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MyMuMu_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MyMuMu() {
  TriggerDictionaryInitialization_MyMuMu_Impl();
}
