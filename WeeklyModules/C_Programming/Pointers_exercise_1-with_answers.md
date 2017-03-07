### Thinking about pointers
Examine the following code examples and describe what happens in each assignment.

```C
char c1 = 'a';
char *char_ptr = &c1; 
```

<!--- 
Declares a variable c1 and initialises with the value 'a'. Declares a pointer to char, char_ptr and assigns the address of c1 
--->

The code above will work perfectly fine. However, the **following code will not work. Why?**

```C
char c1 = 'a';
char *char_ptr;
*char_ptr = &c1; 
```

<!---
```C
char c1 = 'a';
char *char_ptr;
char_ptr = &c1; // Remove the dereference operator
```
--->

Rewrite the above code so that it works correcly

The following code also **will not work**. Why? And how would you fix it?

```C
char c1 = 'a';
char *char_ptr = &c1;
char_ptr = 'b';
```

<!---
```C
char c1 = 'a';
char *char_ptr = &c1;
*char_ptr = 'b'; // Add a dereference operator
```
--->