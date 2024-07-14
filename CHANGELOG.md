# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] - 2024-07-14

### Added
- Support for SIMD on i686

### Fixed
- do not enable `x86_64` specific features on other arches

## [0.2.0] - 2024-07-11

### Added
- Runtime CPU Dispatcher
- configure: `--enable-native` combine compilation units into a single build with all supported isets

### Removed
- configure: `--enable-simd` - removed obsolete option

### Fixed
- Makefile: ensure header changes trigger rebuild
- configure: check for false
- configure: check if pthreads required for `sem_post()`
- pass configure-detected LIBS to LDLIBS

## [0.1.2] - 2023-08-26

Darwin support + minor bugfixes.

### Added
- Darwin support (macOS/aarch64)
- define byteorder macros for macOS/Darwin

### Changed
- tests: prefer arc4random* if available for random bytes

### Fixed
- Ensure tests fail when they exit with an unknown error (eg. SIGSEGV)
- configure: check for <sys/random.h>
- only use getrandom() if header and library found
- minor Makefile fixes

## [0.1.1] - 2023-08-16

Minor bugfixes.

### Added
- added arc4random_uniform() support to rq_random_esi()

### Fixed
- install/uninstall targets use ABI version correctly
- fix realclean make target - remove all built files
- rq_random_esi(): bugfix: ensure value is >= min and <= RQ_ESI_MAX when reading from /dev/urandom
- fix time_t overflow on 32 bit systems
- remove `time` from make speedtest target (not a builtin)
- rq_decode: check return value of rq_encode_symbol
- configure: test if clock_gettime() requires -lrt
- fix warning in example program
- fixing compiler warnings

## [0.1.0] - 2023-07-21

### Added
- OTI common & scheme header functions
    - rq_oti()
    - rq_oti_F()
    - rq_oti_T()
    - rq_oti_Z()
    - rq_oti_N()
    - rq_oti_Al()

### Changed
- updated test runner

### Fixed
- include endian / sys/endian as appropriate
- doc: fix manpage example for rq_decode()
- src/Makefile.in: Directly symlink .so ABI version
    This fixes cross-building on incompatible architectures (e.g. building
    arm packages on x86_64), where ldconfig refuses to create the symlink
    when it does not recognize the non-native architecture.
- only install necessary files from include/
- src/Makefile.in: Use LDFLAGS when building .so file
- if libsodium is enabled, compile liblcrq with an explicit dependency

## [0.0.1] - 2022-07-16

### Added
- bit shifting macros rq_pid2sbn(3), rq_pid2esi(3) rq_pidsetsbn(3) rq_pidsetesi(3)

### Fixed
- CID 274791: Out-of-bounds access
- CID 274784: Unintended sign extension
- CID 274785: Unintended sign extension
- CID 274790: Unintended sign extension
- CID 274793: Unintended sign extension
- CID 274794: Unintended sign extension
- CID 274789: Out-of-bounds read
- CID 274792: Resource leak
- fix big endian bug in rq_symbol(3)
- update examples/ to use renamed API calls

## [0.0.0.0] - 2022-07-14

### Added
- install note for *BSD

### Changed
- remove RFC 6330 from docs due to "non-free" licence to make packaging easier.

## [0.0.0] - 2022-07-12 Initial release

C library RFC6330 RaptorQ Implementation for Librecast
