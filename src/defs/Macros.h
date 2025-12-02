#pragma once
#include <vector>

#define SET_BRANCH_VEC(CPP_TYPE, NAME)         \
    std::vector<CPP_TYPE>* NAME = nullptr;     \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

// Convenience wrappers using ROOT typedefs
#define SET_BRANCH_VFLOAT(NAME)   SET_BRANCH_VEC(Float_t,   NAME)
#define SET_BRANCH_VINT(NAME)     SET_BRANCH_VEC(Int_t,     NAME)
#define SET_BRANCH_VUSHORT(NAME)  SET_BRANCH_VEC(UShort_t,  NAME)
#define SET_BRANCH_VULONG64(NAME) SET_BRANCH_VEC(ULong64_t, NAME)
#define SET_BRANCH_VBOOL(NAME)    SET_BRANCH_VEC(Bool_t,    NAME)
#define SET_BRANCH_VUCHAR(NAME)   SET_BRANCH_VEC(UChar_t,   NAME)

#define SET_BRANCH_ULONG64(NAME) \
    ULong64_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UINT(NAME) \
    UInt_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_INT_ARRAY(NAME) \
    Int_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_FLOAT(NAME) \
    Float_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_FLOAT_ARRAY(NAME) \
    Float_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_BOOL_ARRAY(NAME) \
    Bool_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_BOOL(NAME) \
    Bool_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UCHAR_ARRAY(NAME) \
    UChar_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_USHORT_ARRAY(NAME) \
    UShort_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")
