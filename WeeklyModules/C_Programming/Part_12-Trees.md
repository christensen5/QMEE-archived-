# 12. Trees

Binary trees are an essential part of learning to program in almost any language because they are powerful tools for storing and retrieving data. They are doubly convenient for biologists because they are also useful for representing networks such as phylogenetic trees.

## A review of phylogenetic trees as 'data formats'

### Newick style

The Newick standard for phylogenetic encoding is probably the oldest and most popular format still in general use. It uses a simple system of brackets, commas, and tip names to encode the tree.

Levels of the hierarchy are determined by parentheses, commas separate branches stemming from a node. The standard dictates that tip names are written out in full. However, some encoding conventions allow you to use numeric signifiers. A terminal semicolon indicates the end of the valid Newick string:

`((A,B),(C,D));`

A Newick tree can be preceded by either a `[&R]` or `[&U]` token that indicates whether or not the tree is to be taken as rooted or unrooted. 

### PhyloXML

This is a new and extremely useful standard for encoding phylogenetic information. XML encodings are hierarchical and thus naturally suited to the storage of phylogenetic information. The standards are [documented here](http://www.phyloxml.org/). An example snippet of phyloXML:

```XML
<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns="http://www.phyloxml.org" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd">
    <phylogeny rooted="true">
        <name>Alcohol dehydrogenases</name>
        <description>contains examples of commonly used elements</description>
        <clade>
            <events>
                <speciations>1</speciations>
            </events>
            <clade>
                <taxonomy>
                    <id provider="ncbi">6645</id>
                    <scientific_name>Octopus vulgaris</scientific_name>
                </taxonomy>
                <sequence>
                    <accession source="UniProtKB">P81431</accession>
                    <name>Alcohol dehydrogenase class-3</name>
                </sequence>
            </clade>
            <clade>
                <confidence type="bootstrap">100</confidence>
                <events>
                    <speciations>1</speciations>
                </events>
                <clade>
                    <taxonomy>
                        <id provider="ncbi">1423</id>
                        <scientific_name>Bacillus subtilis</scientific_name>
                    </taxonomy>
                    <sequence>
                        <accession source="UniProtKB">P71017</accession>
                        <name>Alcohol dehydrogenase</name>
                    </sequence>
                </clade>
                <clade>
                    <taxonomy>
                        <id provider="ncbi">562</id>
                        <scientific_name>Escherichia coli</scientific_name>
                    </taxonomy>
                    <sequence>
                        <accession source="UniProtKB">Q46856</accession>
                        <name>Alcohol dehydrogenase</name>
                    </sequence>
                </clade>
            </clade>
        </clade>
    </phylogeny>
</phyloxml>
```
From [here](http://www.phyloxml.org/examples_syntax/phyloxml_syntax_example_1.html)

These formats are useful for storage. They are not, however, good for doing calculations on directly. For that, we can use data structures in C to create trees using areas of memory.

## Trees and nodes in memory

Trees can be created in memory the way linked lists were in a previous section. Thus trees in memory are constructed as a linked set of structures connected by pointers.

![](https://bytebucket.org/mhasoba/silbiocompmasterepo/raw/aa5dcef39a0ef2bab9f48eaf2e54e2d898cea52f/WeeklyModules/C_Programming/images/_ptr_tree.png)

### Defining nodes
Trees are created in memory using structs and pointers. We can start with the following definition

```C
struct node_st {
    struct node_st *left_desc;
    struct node_st *right_desc;
    struct node_st *parent;
};
```

In my programs, because the node struct gets used a lot, I simplify this using a `typedef` statement

```C
typedef struct node_st node_t;

typedef struct node_st {
    node_t *left_desc;
    node_t *right_desc;
    node_t *parent;
} node_t;
```

The typedef statement therefore just allows me to abbreviate `struct node_st` to just `node_t`.

So, I have now created a new variable type called `node_t`. You can expand this struct with more internal variables as and when you need.

For instance, we might need a variable to indicate whether the node is a tip.

```C
typedef struct node_st {
    node_t *left_desc;
    node_t *right_desc;
    node_t *parent;
    int tip;
} node_t;
```

### Defining trees

Trees are based on nodes. Thus, we can assemble a tree from a set of node structures. Thus, at the core of a tree would be a block of node structures:

```C
node_t *treenodes;
```

We could then use `malloc` or `calloc` to give us exactly as many nodes as we needed. However, this doesn't give us the ability to create and destroy nodes as we need them. An alternative would be to create a block of pointers to node structures, and allocate the node memory individually. We could, in this instance, create and destroy nodes as we need them.

```C
node_t **treenodes;
```

The overall tree structure should include all the parameters we need to safely work with the tree

```C
typedef struct tree_st {
	int n_spp;			// The number of tips
	int n_nodes;		// The number of internal nodes
	/* Extend this structure by including your new variables here
	*/
	node_t **treenodes;	// For the array of nodes
	node_t *start;		// For the root or start node
} tree_t;
```

## Traversing the tree

Once we have designed a tree structure, we can design a function that would traverse a tree and which could be used to apply **preorder** or **postorder** functions to the tree:

```C
void basic_tree_traversal(node_t* n)
{
    if (n->tip) {
        return;
    }

    basic_tree_traversal(n->left_desc);
    basic_tree_traversal(n->right_desc);
}
```

We call such a function on the pointer to the root node of the tree. **A traversal of this kind on a completely bifurcating tree is guaranteed to visit every tip and every node exactly once.**


# Practical problem

Write a program that reads a Newick-style tree encoding and builds a tree in memory.

1. Design the data types you need to use (we can start with what we have in this file)
2. Define these
3. Think about how you will 'read' the Newick tree and how that will translate to a tree in memory. What needs to happen first?
4. Start with some simple functions. Perhaps something that prints a Newick tree 'clade by clade'. 
5. Then work up from there.