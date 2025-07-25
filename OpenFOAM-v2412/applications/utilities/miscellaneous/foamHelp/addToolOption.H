argList::addArgument("tool");

argList::addNote("Valid <tool> options include:");
for (const word& tool : helpType::dictionaryConstructorTablePtr_->sortedToc())
{
    argList::addNote("    " + tool);
}

argList::addNote
(
    "\nNOTE the <tool> must actually appear *before* any options"
);
