    // Construct optimisation manager
    autoPtr<optimisationManager> optManagerPtr
    (
        optimisationManager::New(mesh)
    );
    optimisationManager& om = optManagerPtr();
