import { Lightbulb, Filter, Target, ListOrdered } from 'lucide-react';

export function Uses() {
  const pipelines = [
    {
      icon: Lightbulb,
      title: 'Generative Screening Pipeline',
      purpose: 'Invents entirely new molecules similar to provided seed compounds using deep generative AI.',
      howItWorks: 'Upload a CSV of promising molecules as seeds; the pipeline learns chemical patterns and generates novel structures. The predicted properties for each AI-generated molecule, such as activity and safety, are included in the output.',
      why: 'Explores chemical space humans might never consider; ideal for early-stage hit discovery; generates candidates for further screening.',
      gradient: 'from-purple-500 to-pink-500'
    },
    {
      icon: Filter,
      title: 'Virtual Screening Pipeline',
      purpose: 'Filters massive libraries (thousands+) to identify the most promising drug candidates.',
      howItWorks: 'Upload a large CSV of candidates; AI ranks molecules based on predicted performance, such as activity against a target. Download a filtered set of high-scoring hits to test first.',
      why: 'Focuses resources on the best candidates, dramatically reducing costly experimental trial and error.',
      gradient: 'from-purple-500 to-blue-500'
    },
    {
      icon: Target,
      title: 'Predict Properties Pipeline',
      purpose: 'Estimates properties (activity, solubility, toxicity, drug-likeness) for submitted molecules.',
      howItWorks: 'Upload any set of molecules; the pipeline outputs predictions as new columns—early feedback on likelihood of success, safety, and physical properties.',
      why: 'Allows rapid hypothesis testing; helps decide whether to move forward or revisit designs before lab work starts.',
      gradient: 'from-purple-500 to-indigo-500'
    },
    {
      icon: ListOrdered,
      title: 'Prioritize Candidates Pipeline',
      purpose: 'Ranks candidate molecules so you know which ones to synthesize or test first, based on predicted properties.',
      howItWorks: 'Upload a CSV of candidates with predicted properties; pipeline sorts and prioritizes based on your chosen criteria (activity, safety, drug-likeness). Output includes a priority ranking column.',
      why: 'Ensures best use of lab time and resources; maximizes chance of success by focusing on the most promising compounds.',
      gradient: 'from-purple-500 to-violet-500'
    }
  ];

  return (
    <div className="max-w-7xl mx-auto space-y-12 animate-in fade-in duration-700">
      {/* Header */}
      <div className="text-center space-y-4">
        <h1 className="text-5xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent">
          Pipeline Uses
        </h1>
        <p className="text-xl text-gray-400">
          Four powerful AI pipelines to accelerate your drug discovery workflow
        </p>
      </div>

      {/* Pipelines */}
      <div className="space-y-8">
        {pipelines.map((pipeline, index) => {
          const Icon = pipeline.icon;
          return (
            <div
              key={index}
              className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-2xl p-8 backdrop-blur-sm hover:border-purple-400/50 transition-all hover:shadow-xl hover:shadow-purple-500/20"
            >
              {/* Title Section */}
              <div className="flex items-start gap-6 mb-6">
                <div className={`p-4 bg-gradient-to-br ${pipeline.gradient} rounded-xl shadow-lg`}>
                  <Icon className="w-8 h-8 text-white" />
                </div>
                <div className="flex-1">
                  <h2 className="text-3xl text-white mb-2">{pipeline.title}</h2>
                  <div className="w-24 h-1 bg-gradient-to-r from-purple-500 to-purple-700 rounded-full"></div>
                </div>
              </div>

              {/* Content Grid */}
              <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                <div className="space-y-2">
                  <div className="inline-block px-3 py-1 bg-purple-500/20 rounded-full text-sm text-purple-300 mb-2">
                    Purpose
                  </div>
                  <p className="text-gray-300">{pipeline.purpose}</p>
                </div>

                <div className="space-y-2">
                  <div className="inline-block px-3 py-1 bg-blue-500/20 rounded-full text-sm text-blue-300 mb-2">
                    How it Works
                  </div>
                  <p className="text-gray-300">{pipeline.howItWorks}</p>
                </div>

                <div className="space-y-2">
                  <div className="inline-block px-3 py-1 bg-green-500/20 rounded-full text-sm text-green-300 mb-2">
                    Why
                  </div>
                  <p className="text-gray-300">{pipeline.why}</p>
                </div>
              </div>
            </div>
          );
        })}
      </div>

      {/* Call to Action */}
      <div className="bg-gradient-to-r from-purple-900/40 to-purple-800/40 border border-purple-500/40 rounded-2xl p-10 text-center">
        <h3 className="text-3xl text-white mb-4">Ready to accelerate your research?</h3>
        <p className="text-gray-300 mb-6">
          Start using our AI pipelines today and discover molecules faster than ever before
        </p>
        <button className="px-8 py-3 rounded-full bg-gradient-to-r from-purple-500 to-purple-700 text-white hover:shadow-lg hover:shadow-purple-500/50 transition-all">
          Get Started
        </button>
      </div>
    </div>
  );
}
