# summary
This is the spot for code that can be considered a "general" numerical types and methods.

I find it useful to separate such methods based on their numerical domain which is defined by 2 components:
1. The data type (e.g. double, float, int, ...)
2. The spatial dimensionality (e.g. any, 1, 2, ...)

For data types, I use a short identifier
* `D` = double = double-precision floating point numbers
* `F` = float  = single-precision floating point numbers
* `I` = int    = 32-bit integers
* `L` = long   = 64-bit integers

I create short namespace names (with a corresponding folder structure) as follows.
* `DX` -> Types and methods related to 64-bit floating-point numbers of arbitrary spatial dimensionality.
* `D1` -> Types and methods related to 64-bit floating-point numbers in 1 spatial dimension.
* `D2` -> Types and methods related to 64-bit floating-point numbers in 2 spatial dimensions.
* `D2` -> Types and methods related to 32-bit floating-point numbers in 2 spatial dimensions.
* `D3` -> Types and methods related to 32-bit floating-point numbers in 3 spatial dimensions.

So for example, I could create 3 different data types, each of which represents a spatial coordinate.  I could
name them all `Point`and use the namespace (folder) to indicate the difference.
* `D2::Point`
* `F2::Point`
* `F3::Point`

It is helpful for a human reader to see inline types written as something like `D2::Point` rather than `Point`
because it resolves ambiguity.  For this reason, all inline references to these types should include the short
namespace.
