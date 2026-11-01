# CGV project rules

1. **Purpose**: a fast, light, secure BAM/CRAM viewer for medical and human genetics, replacing
   IGV; online, desktop and company-server versions. Clinical scenarios and QC are out of scope for
   now.
2. **Priorities**: speed and security.
3. **Security**: our own code, with no dependencies that bring risk. Third-party libraries are
   welcome when they are safe and make development simpler — taken as a reviewed copy at a fixed
   version, never installed from a package registry. Our own code must be free of threats too:
   replacing old code with risky dependencies by new code that is still unsafe is not acceptable.
   Details and checks: SECURITY.md.
4. **Work in stages**, each checked by the user before the next one.
5. **Interface**: as compact and information-dense as possible — the screen goes to data. Short
   labels, small controls, no decorative space. (Confirmed 2026-11-01: so far the interface is on
   the right track.)
6. **Like IGV where IGV is good**: users come from IGV, so behavior and terms follow it unless
   there is a reason to differ.
7. **Language**: interface and documents in English by default; Russian as an option (interface
   switch, `docs/*.ru.md`).
8. **Documentation**: the user guide (`docs/USER_GUIDE.md`, `docs/USER_GUIDE.ru.md`) describes
   every function as it works and looks now, and is updated with every change. Explanations belong
   in the guide; in the interface only as tooltips, where possible and without overloading it (the
   browser's own tooltip, shown after the pointer rests, never on every pass; none over the reads).
   Anything CGV shows that is not a common standard (e.g. mQV) is explained in the guide.
   Future work goes to `docs/TODO.md`.
