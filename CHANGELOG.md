# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added
- bit shifting macros rq_pid2sbn(3), rq_pid2esi(3) rq_pidsetsbn(3) rq_pidsetesi(3)

### Fixed
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
