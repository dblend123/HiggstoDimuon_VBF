#pragma once

#define SET_BRANCH_ULONG64(NAME) \
    ULong64_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UINT(NAME) \
    UInt_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_FLOAT(NAME) \
    Float_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

    /// TODO
    /// customize further
#define SET_BRANCH_INT_ARRAY(NAME) \
    Int_t NAME[100]; \
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

#define SET_BRANCH_UCHAR_ARRAY(NAME) \
    UChar_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_USHORT_ARRAY(NAME) \
    UShort_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")
