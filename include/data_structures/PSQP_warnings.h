#ifndef PSQP_WARNINGS_H
#define PSQP_WARNINGS_H

#if defined(__clang__)
#define PSQP_DIAG_PUSH() _Pragma("clang diagnostic push")
#define PSQP_DIAG_POP() _Pragma("clang diagnostic pop")

#define PSQP_DIAG_IGNORE_CONVERSION()                                               \
    _Pragma("clang diagnostic ignored \"-Wconversion\"")

#define PSQP_DIAG_IGNORE_SIGN_CONVERSION()                                          \
    _Pragma("clang diagnostic ignored \"-Wsign-conversion\"")

#elif defined(__GNUC__)
#define PSQP_DIAG_PUSH() _Pragma("GCC diagnostic push")
#define PSQP_DIAG_POP() _Pragma("GCC diagnostic pop")

#define PSQP_DIAG_IGNORE_CONVERSION()                                               \
    _Pragma("GCC diagnostic ignored \"-Wconversion\"")

#define PSQP_DIAG_IGNORE_SIGN_CONVERSION()                                          \
    _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"")

#elif defined(_MSC_VER)
#define PSQP_DIAG_PUSH() __pragma(warning(push))
#define PSQP_DIAG_POP() __pragma(warning(pop))

/* Closest MSVC equivalents */
#define PSQP_DIAG_IGNORE_CONVERSION() __pragma(warning(disable : 4244 4267))

#define PSQP_DIAG_IGNORE_SIGN_CONVERSION() __pragma(warning(disable : 4245))

#else
#define PSQP_DIAG_PUSH()
#define PSQP_DIAG_POP()
#define PSQP_DIAG_IGNORE_CONVERSION()
#define PSQP_DIAG_IGNORE_SIGN_CONVERSION()
#endif

#endif /* PSQP_WARNINGS_H */
