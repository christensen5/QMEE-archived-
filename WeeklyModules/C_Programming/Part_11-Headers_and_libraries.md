# 11. Header files, programs with multiple files, and preprocessor statements

We will not to into great detail on these topics here. You should read up on these topics before attempting to create any larger programs in C. However, for relatively simple tasks, like outsourcing some computationally 'expensive' inner-loop routines to a set of C functions, you won't need to worry about these. However, you should be aware of them, what they are and basically what they do. Header file inclusion is important in communicating between modules written in different languages. Secondly, the preferred interface for calling C code from R makes extensive use of macros that, without introduction, will look really alien to you.

## Header files

Header files are indicated with a `.h` extension (as we've seen) and contain various definitions, function prototypes, variable and data type definitions that are required throughout a program or library. We've been including multiple header files already. Most of these are standard C libraries that live on pre-defined 'search paths' on your machine. When including files on the pre-defined search paths, the file name is included in angled brackets:

```C
#include <stdbool.h>
```

However, it is customary to include many definitions and function prototypes in a separate header file. This allows you to more easily work with programs involving multiple files. Therefore, your local directory for source files might have its own header files that you created. To instruct the pre-processor that you expect those files to be living in the same directory as the source code, you simply wrap the name in double quotes:

```C
#include "funclibr.h"
```

## Compiling multiple files

It is common to require multiple files for a program. This helps organise logically connected 'chunks' of source code and is extremely useful when more than one person is working on the program's development. If your source code consists of multiple files, then you compile these simply by including multiple input files in the 

```bash
gcc main.c funclibr.c
```

However, in order for functions in `main.c` to be able to use functions in `funclib.c`, the `main.c` file will need to have the prototypes of those functions. Therefore, we can ensure that all of a file's function prototypes live in its header file that gets included as appropriate.

## Define and Macros

We won't cover these in detail, but mainly mention their existence and a bit about how they work. Macros and `define` statements allow us to instruct the preprocessor, prior to compiling, that there are some 'shorthands' we have defined for other values. Like other preprocessor directives (for instance the `include` keyword), preporcessor definitions start with the `#` symbol. Preprocessor directives extend to the end of the line. Thus if you need to embed newlines within your macro definitions, you need to 'escape' them with a backslash.

### the `define` keyword

The `define` keyword lets us substitute values for a symbolic representation.  For instance, we might have some numeric values that we will use frequently but don't want to type out all the time:

```C
#define BASE_NAT_LOG 2.71828
```

Now, anywhere we write `BASE_NAT_LOG` in our code will expand to `2.71828`

It is not unheard of for programmers to use preprocessor definitions to silently 'sweep up' common typos:

`#define pritnf printf`

This practice is not recommended.

Typically, you would use a preprocessor macro to improve the portability and extendabiliy of your code. Instead of using 'magic numbers', you will save yourself a lot of headaches by using a preprocessor definition in your code. That way, if you need to change a constant value, you don't have to search your entire source code to change all instances. You can just change the definition.

Note also, that it is common to create preprocessor definitions in all-caps to distinguish them from other variables.

### Preprocessor macros

Preprocessor macros allow us to substitute large pieces of code for just a simple short word. These preprocessor macros can even take arguments. Try the following program and see if you can figure out how this #define statment works (and why I needed to wrap everything in braces).

```C
#include <stdio.h>

#define MIN(x, y, result) {if (x > y){ \
					result = y; \
				} \
				else {\
					result = x; \
				}}

int main (void)
{

	int a = 2;
	int b = 3;
	int c = 0;

	MIN(a, b, c);

	printf("Minimum: %i\n", c);

	return 0;
}
```

# Exercises

### 1- Safe arrays
Create a library and header file for 'safe arrays'.

### 2- Advanced exercise: 
Most modern machines only allow up to 64-bit integer widths. However, one may wish to exploit bitwise operations that require larger sets of bits. Support for 128-bit integers is a bit ambiguous, and possibly only available with some compilers. Create a library for doing bitwise operations on bit sets of arbitrary size. How would you design such a library? What features of C would you use?