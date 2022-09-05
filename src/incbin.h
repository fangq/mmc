#pragma once

#define GET_MACRO(_1,_2,_3,_4,NAME,...) NAME

#ifdef __EMSCRIPTEN__
#define INCBIN_4(type, symbol, length, file) \
    __asm( \
           ".section .data, \"\", @\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           ".size " #symbol ", .-" #symbol "\n" \
           #length ":\n" \
           ".long .-" #symbol "\n" \
           ".size " #length ", 4\n"); \
    extern const type symbol[] __asm("" #symbol); \
    extern const unsigned int length  __asm("" #length);

#define INCBIN_2(symbol, file) INCBIN(unsigned char, symbol, symbol ## Length, file)

#define INCBIN(...) GET_MACRO(__VA_ARGS__, INCBIN_4, NULL, INCBIN_2)(__VA_ARGS__)

#define INCTXT_3(symbol, length, file) \
    __asm( \
           ".section .data, \"\", @\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           ".byte 0\n" \
           ".size " #symbol ", .-" #symbol "\n" \
           #length ":\n" \
           ".long .-" #symbol "\n" \
           ".size " #length ", 4\n"); \
    extern const char symbol[] __asm("" #symbol); \
    extern const unsigned int length  __asm("" #length);

#define INCTXT_2(symbol, file) \
    __asm( \
           ".section .data, \"\", @\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           ".byte 0\n" \
           ".size " #symbol ", .-" #symbol "\n"); \
    extern const char symbol[] __asm("" #symbol);

#define INCTXT(...) GET_MACRO(__VA_ARGS__, NULL, INCTXT_3, INCTXT_2)(__VA_ARGS__)

#else

#define INCBIN_4(type, symbol, length, file) \
    __asm( \
           ".section .data\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           #length ":\n" \
           ".long .-" #symbol "\n"); \
    extern const type symbol[] __asm("" #symbol); \
    extern const unsigned int length  __asm("" #length);

#define INCBIN_2(symbol, file) INCBIN(unsigned char, symbol, symbol ## Length, file)

#define INCBIN(...) GET_MACRO(__VA_ARGS__, INCBIN_4, NULL, INCBIN_2)(__VA_ARGS__)

#define INCTXT_3(symbol, length, file) \
    __asm( \
           ".section .data\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           ".byte 0\n" \
           #length ":\n" \
           ".long .-" #symbol "\n"); \
    extern const char symbol[] __asm("" #symbol); \
    extern const unsigned int length  __asm("" #length);

#define INCTXT_2(symbol, file) \
    __asm( \
           ".section .data\n" \
           #symbol ":\n" \
           ".incbin " #file "\n" \
           ".byte 0\n"); \
    extern const char symbol[] __asm("" #symbol);

#define INCTXT(...) GET_MACRO(__VA_ARGS__, NULL, INCTXT_3, INCTXT_2)(__VA_ARGS__)

#endif