# 12. Trees

Binary trees are an essential part of learning to program in almost any language because they are powerful tools for storing and retrieving data. They are doubly convenient for biologists because they are also useful for representing networks such as phylogenetic trees.

## A review of phylogenetic trees as 'data formats'

### Newick style



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
