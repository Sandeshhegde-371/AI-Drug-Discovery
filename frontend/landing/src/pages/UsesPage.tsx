import { Lightbulb, Filter, Target, ListOrdered, Home, Brain, ArrowRight } from 'lucide-react';
import { Link } from 'react-router-dom';

export function UsesPage() {
  const pipelines = [
    {
      icon: Lightbulb,
      title: 'Generative Screening Pipeline',
      purpose: 'Invents entirely new molecules similar to provided seed compounds using deep generative AI.',
      howItWorks: 'Upload a CSV of promising molecules as seeds; the pipeline learns chemical patterns and generates novel structures. The predicted properties for each AI-generated molecule, such as activity and safety, are included in the output.',
      why: 'Explores chemical space humans might never consider; ideal for early-stage hit discovery; generates candidates for further screening.',
      gradient: 'from-purple-500 to-pink-500',
      bgGradient: 'from-purple-950/40 to-pink-950/20',
      borderColor: 'border-purple-500/30',
      hoverBorder: 'hover:border-purple-400/50',
      badgeColors: {
        purpose: 'bg-purple-500/20 text-purple-300',
        howItWorks: 'bg-pink-500/20 text-pink-300',
        why: 'bg-purple-400/20 text-purple-200'
      }
    },
    {
      icon: Filter,
      title: 'Virtual Screening Pipeline',
      purpose: 'Filters massive libraries (thousands+) to identify the most promising drug candidates.',
      howItWorks: 'Upload a large CSV of candidates; AI ranks molecules based on predicted performance, such as activity against a target. Download a filtered set of high-scoring hits to test first.',
      why: 'Focuses resources on the best candidates, dramatically reducing costly experimental trial and error.',
      gradient: 'from-blue-500 to-cyan-500',
      bgGradient: 'from-blue-950/40 to-cyan-950/20',
      borderColor: 'border-blue-500/30',
      hoverBorder: 'hover:border-blue-400/50',
      badgeColors: {
        purpose: 'bg-blue-500/20 text-blue-300',
        howItWorks: 'bg-cyan-500/20 text-cyan-300',
        why: 'bg-blue-400/20 text-blue-200'
      }
    },
    {
      icon: Target,
      title: 'Predict Properties Pipeline',
      purpose: 'Estimates properties (activity, solubility, toxicity, drug-likeness) for submitted molecules.',
      howItWorks: 'Upload any set of molecules; the pipeline outputs predictions as new columns—early feedback on likelihood of success, safety, and physical properties.',
      why: 'Allows rapid hypothesis testing; helps decide whether to move forward or revisit designs before lab work starts.',
      gradient: 'from-indigo-500 to-purple-500',
      bgGradient: 'from-indigo-950/40 to-purple-950/20',
      borderColor: 'border-indigo-500/30',
      hoverBorder: 'hover:border-indigo-400/50',
      badgeColors: {
        purpose: 'bg-indigo-500/20 text-indigo-300',
        howItWorks: 'bg-purple-500/20 text-purple-300',
        why: 'bg-indigo-400/20 text-indigo-200'
      }
    },
    {
      icon: ListOrdered,
      title: 'Prioritize Candidates Pipeline',
      purpose: 'Ranks candidate molecules so you know which ones to synthesize or test first, based on predicted properties.',
      howItWorks: 'Upload a CSV of candidates with predicted properties; pipeline sorts and prioritizes based on your chosen criteria (activity, safety, drug-likeness). Output includes a priority ranking column.',
      why: 'Ensures best use of lab time and resources; maximizes chance of success by focusing on the most promising compounds.',
      gradient: 'from-violet-500 to-fuchsia-500',
      bgGradient: 'from-violet-950/40 to-fuchsia-950/20',
      borderColor: 'border-violet-500/30',
      hoverBorder: 'hover:border-violet-400/50',
      badgeColors: {
        purpose: 'bg-violet-500/20 text-violet-300',
        howItWorks: 'bg-fuchsia-500/20 text-fuchsia-300',
        why: 'bg-violet-400/20 text-violet-200'
      }
    }
  ];

  return (
    <div className="min-h-screen bg-[#0a0a0f] text-white">
      {/* Gradient Backgrounds */}
      <div className="fixed inset-0 overflow-hidden pointer-events-none">
        <div className="absolute top-0 right-0 w-[600px] h-[600px] bg-blue-600/20 rounded-full blur-[120px]" />
        <div className="absolute bottom-0 left-0 w-[500px] h-[500px] bg-purple-600/15 rounded-full blur-[100px]" />
      </div>

      <div className="relative z-10">
        {/* Header */}
        <header className="px-8 py-6 flex items-center justify-between backdrop-blur-sm bg-white/5 border-b border-white/10 sticky top-0 z-50">
          <Link to="/" className="flex items-center gap-3 hover:opacity-80 transition-opacity group">
            <div className="p-2 bg-blue-500/20 rounded-lg group-hover:bg-blue-500/30 transition-colors">
              <Home className="w-5 h-5 text-blue-400" />
            </div>
            <span className="text-white/80">Back to Home</span>
          </Link>
          
          <div className="flex items-center gap-3">
            <div className="p-2 bg-gradient-to-br from-blue-500 to-blue-700 rounded-lg">
              <Brain className="w-6 h-6 text-white" />
            </div>
            <div className="flex flex-col">
              <span className="text-xl bg-gradient-to-r from-blue-400 to-blue-600 bg-clip-text text-transparent">
                AI-DRIVEN DDD
              </span>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            <Link to="/overview" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Overview
            </Link>
            <Link to="/sample" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Sample
            </Link>
          </div>
        </header>

        {/* Content */}
        <main className="container mx-auto px-8 py-16">
          <div className="max-w-7xl mx-auto space-y-16">
            {/* Header */}
            <div className="text-center space-y-6">
              <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-blue-500/10 border border-blue-500/20 mb-4">
                <Filter className="w-4 h-4 text-blue-400" />
                <span className="text-sm text-blue-300">AI Pipelines</span>
              </div>
              
              <h1 className="text-6xl leading-tight">
                <span className="block text-white mb-2">Powerful</span>
                <span className="bg-gradient-to-r from-blue-400 via-purple-500 to-purple-600 bg-clip-text text-transparent">
                  AI Pipelines
                </span>
              </h1>
              
              <p className="text-xl text-gray-400 max-w-3xl mx-auto">
                Four specialized pipelines to accelerate your drug discovery workflow
              </p>
            </div>

            {/* Pipelines */}
            <div className="space-y-10">
              {pipelines.map((pipeline, index) => {
                const Icon = pipeline.icon;
                return (
                  <div
                    key={index}
                    className={`relative bg-gradient-to-br ${pipeline.bgGradient} border ${pipeline.borderColor} rounded-3xl p-10 backdrop-blur-sm ${pipeline.hoverBorder} transition-all hover:shadow-2xl overflow-hidden group`}
                  >
                    {/* Decorative glow */}
                    <div className="absolute -top-24 -right-24 w-48 h-48 bg-gradient-to-br opacity-0 group-hover:opacity-20 transition-opacity rounded-full blur-3xl" 
                         style={{background: `linear-gradient(to bottom right, var(--tw-gradient-stops))`}} />
                    
                    <div className="relative z-10">
                      {/* Title Section */}
                      <div className="flex items-start gap-6 mb-8">
                        <div className={`p-5 bg-gradient-to-br ${pipeline.gradient} rounded-2xl shadow-lg group-hover:scale-110 transition-transform`}>
                          <Icon className="w-10 h-10 text-white" />
                        </div>
                        <div className="flex-1">
                          <h2 className="text-4xl text-white mb-3">{pipeline.title}</h2>
                          <div className={`w-32 h-1.5 bg-gradient-to-r ${pipeline.gradient} rounded-full`}></div>
                        </div>
                        <div className="hidden md:block px-4 py-2 bg-white/5 rounded-full border border-white/10">
                          <span className="text-sm text-gray-400">Pipeline {index + 1}</span>
                        </div>
                      </div>

                      {/* Content Grid */}
                      <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
                        <div className="space-y-3">
                          <div className={`inline-block px-4 py-1.5 ${pipeline.badgeColors.purpose} rounded-full text-sm mb-2`}>
                            Purpose
                          </div>
                          <p className="text-gray-300 leading-relaxed">{pipeline.purpose}</p>
                        </div>

                        <div className="space-y-3">
                          <div className={`inline-block px-4 py-1.5 ${pipeline.badgeColors.howItWorks} rounded-full text-sm mb-2`}>
                            How it Works
                          </div>
                          <p className="text-gray-300 leading-relaxed">{pipeline.howItWorks}</p>
                        </div>

                        <div className="space-y-3">
                          <div className={`inline-block px-4 py-1.5 ${pipeline.badgeColors.why} rounded-full text-sm mb-2`}>
                            Why Use It
                          </div>
                          <p className="text-gray-300 leading-relaxed">{pipeline.why}</p>
                        </div>
                      </div>
                    </div>
                  </div>
                );
              })}
            </div>

            {/* Call to Action */}
            <div className="bg-gradient-to-r from-purple-900/40 to-blue-800/40 border border-purple-500/40 rounded-3xl p-12 text-center">
              <h3 className="text-4xl text-white mb-4">Ready to accelerate your research?</h3>
              <p className="text-gray-300 mb-8 text-lg max-w-2xl mx-auto">
                Start using our AI pipelines today and discover molecules faster than ever before
              </p>
              <Link 
                to="/sample"
                className="inline-flex items-center gap-3 px-8 py-4 rounded-full bg-gradient-to-r from-purple-500 to-blue-600 text-white hover:shadow-lg hover:shadow-purple-500/50 transition-all hover:scale-105"
              >
                <span>View Sample Data</span>
                <ArrowRight className="w-5 h-5" />
              </Link>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
}
