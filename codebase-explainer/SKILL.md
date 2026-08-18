---
name: codebase-explainer
description: Explain how code works — whole codebases, specific features or flows through a codebase, or individual functions/classes/files. Use this skill whenever the user asks how something in a codebase works, what a file or function does, how data flows through code, how a feature is implemented, or asks for an architectural tour or orientation in an unfamiliar codebase. Trigger even when they don't use the word "explain" — questions like "what does this do", "how is X implemented", "walk me through this", "what's going on in here", or "I'm new to this repo, where do I start" all count.
---

# Codebase Explainer

A skill for explaining code to a strong engineer who wants conceptual understanding, not narration.

## The core idea

The user wants to understand how code *works conceptually*, not have it read back to them. The most common failure mode is going line-by-line or block-by-block saying "then this happens, then this happens" — a zoomed-in narration of the source rather than an explanation of what the code is *doing* and *why*. Avoid this above all else.

Lead with the concept, then ground it in the code. The concept is the data flow, the responsibility of each part, the relationships between pieces. The code is where those concepts live — useful as reference, but not the explanation itself.

## What kind of question is this

Infer from the question which type of explanation is wanted — don't ask. The three types are different tasks with different shapes of answer:

- **Whole codebase tour** — orientation in an unfamiliar repo. "What is this?", "walk me through this codebase", "I'm new to this repo".
- **Feature or flow trace** — how a particular thing works end-to-end. "How does authentication work?", "how does X get from the API to the database?", "what happens when a user does Y?".
- **Zoom-in on a unit** — what a specific function, class, or file does. "What does this function do?", "explain this class".

The question almost always makes it obvious which one it is. In rare genuinely-ambiguous cases, pick the most likely interpretation and answer it — the user will redirect if needed. Asking a clarifying question first is more annoying than guessing wrong.

## Always do this first

Investigate before answering. The biggest source of wrong explanations is answering from the code that happens to be in view without checking how it connects to the rest. Before composing the answer:

- For a **whole codebase tour**: skim the top-level structure, the README if there is one, the main entry points, and the central model/type definitions. Get a sense of what kind of project it is (library, web app, CLI, pipeline, etc.).
- For a **feature/flow trace**: find all the relevant pieces. Where does the flow start? Where does it end? What gets called in between? Use grep/search to find call sites and definitions rather than guessing.
- For a **zoom-in on a unit**: read the unit itself, and also find where it's called from. The caller's context often makes the function's purpose obvious in a way the function alone doesn't.

If the investigation is partial — you looked at the models but not the views, or you traced the happy path but not the error handling — say so when you answer. An honest "I didn't look at X" is more useful than a confident summary that turns out to be missing the important bit.

## The code is the source of truth

Base the explanation on what the code actually does, not what the documentation says it does. READMEs, docstrings, and design docs go out of date — they describe what the project was, or what someone meant for it to be, not necessarily what it is now. If the README says the library does X but the code clearly does Y, the answer is Y.

Use documentation only as a supporting signal: to confirm something the code already shows, to clarify intent where the code is ambiguous, or to find pointers to where things live. Never let a doc override what the code is plainly doing. If there's a contradiction worth flagging (e.g. the README's description is materially misleading), mention it — but the explanation itself should describe the code's actual behaviour.

## Shape of the answer

### Whole codebase tour

Cover, roughly in this order:

1. **What it is.** One or two sentences on what the project does. Not "this is a Python codebase" — what does it actually accomplish?
2. **Key concepts and object types/models.** The central nouns of the codebase. In a Django app these are the models; in a parser they might be the AST node types; in a pipeline they're the data structures that flow between stages. Name them and explain what each represents. This is usually the most important part of the tour — once the user knows the central types, a lot of the rest falls into place.
3. **Running processes** (if it's an application that runs as a service or system). Enumerate everything that runs: web server, background workers, schedulers, daemons, cron jobs, sidecars. For each, what it does and what it talks to (database, message broker, cache, external API). If it's a library or CLI rather than a running system, skip this.
4. **Top-level structure.** The shape of the repo at the top — the main directories and what each is for. Not a file-by-file listing, but a sense of how the codebase is organised: is it split by layer (models/views/services), by domain (billing/auth/reports), by component (frontend/backend/worker), monorepo with separate packages, etc. The user should come away knowing what kind of organisation to expect when navigating.

### Feature or flow trace

Trace the data, not the source files. The shape is usually: input arrives somewhere → gets transformed/validated → hits some core logic → produces a result → goes somewhere. Describe that flow in prose, naming the functions/classes/methods that handle each stage as you go.

If a function in the flow calls another function that does meaningful work, explain what that does too — don't leave a black box in the middle of the trace. But don't recurse indefinitely: if function A calls B calls C calls D, summarise the deeper levels rather than tracing every one. The user wants to understand the flow, not see a complete call graph.

### Zoom-in on a unit

Explain what the unit's job is — its responsibility, what it takes in, what it produces, what side effects it has. Then explain how it does that job, conceptually. If the implementation has any non-obvious parts (a clever trick, a workaround, a subtle invariant), call those out.

If the unit is called from elsewhere, mention where and why — that context often clarifies what the unit is *for*. But again, don't recurse: if the caller is itself called by something complex, stop there.

## Things to avoid

**Don't narrate the source.** "First we declare a variable, then we loop over the items, then we call save..." is reading code aloud. The user can read the code themselves. Explain what the code is *doing*, not what each line *is*.

**Don't use analogies.** They tend to be more confusing than helpful and often misrepresent the actual mechanics. Explain the thing directly.

**Don't assume knowledge of other parts of the codebase.** If you reference `UserManager` or `the validation pipeline` or `the standard flow`, briefly say what that is. The user might be looking at this code for the first time and shouldn't have to ask follow-up questions to decode the explanation.

**Don't lead with line numbers and file paths.** Names of types, functions, and concepts are what carry meaning. "The `Component` model owns stereochemistry via its `StereoGroup` relation" is useful; "in `models.py` at line 47 there's a class" is not. Locations are fine as supporting detail, but the explanation should make sense without them.

**Don't pad with generic framework facts.** "Django uses an ORM" or "FastAPI is an async web framework" is rarely what the user wants to hear. Engage with what's actually specific or interesting about *this* codebase.

## Diagrams

Don't use diagrams by default. The execution environment may not render them — a terminal won't show mermaid at all, and many editor agents don't either. A raw block of mermaid syntax dumped into a non-rendering surface is worse than no diagram.

Only use a diagram if the user explicitly asks for one, or if the rendering environment is known to support it. Otherwise, explain in prose.

## Output

Respond conversationally in the current surface. Don't produce a markdown file, document, or other artifact unless the user explicitly asks for one.